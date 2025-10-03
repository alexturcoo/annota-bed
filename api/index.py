# api/index.py
import os, io, math, tempfile, uuid, shutil
from collections import defaultdict, Counter

from flask import Flask, request, jsonify, send_file
import pandas as pd
import requests

app = Flask(__name__)

# --- Try to import pysam only if available (works locally / on external host) ---
try:
    import pysam
    HAS_PYSAM = True
except Exception:
    HAS_PYSAM = False

FLASK_UPSTREAM = os.getenv("FLASK_UPSTREAM")  # e.g., https://your-flask-host.com

# ----------------------- Shared helpers -----------------------
def parse_bed(file_bytes: bytes):
    out = []
    for line in io.BytesIO(file_bytes).read().decode(errors="ignore").splitlines():
        if not line or line.startswith(("track", "browser", "#")):
            continue
        parts = line.split()
        if len(parts) < 3:
            continue
        chrom = parts[0]
        if chrom.startswith("chr"):
            chrom = chrom[3:]
        try:
            start = int(parts[1]); end = int(parts[2])
        except ValueError:
            continue
        if end > start:
            out.append((chrom, start, end))
    return out

def json_error(msg, code=400):
    return jsonify({"error": msg}), code

# ----------------------- Health -----------------------
@app.get("/health")
def health():
    mode = "proxy->upstream" if FLASK_UPSTREAM else ("full-pysam" if HAS_PYSAM else "serverless-demo")
    return {"status": "ok", "mode": mode}

# ----------------------- Annotate -----------------------
@app.post("/annotate")
def annotate():
    # 0) If we have an upstream, forward the exact multipart request and return its response
    if FLASK_UPSTREAM:
        try:
            # Forward multipart form to upstream /annotate
            files = {}
            for key in ("file", "gtf_bgz", "gtf_tbi"):
                f = request.files.get(key)
                if f:
                    files[key] = (f.filename, f.stream, f.mimetype or "application/octet-stream")
            data = dict(request.form)  # ambiguities, etc.

            resp = requests.post(f"{FLASK_UPSTREAM.rstrip('/')}/annotate", files=files, data=data, timeout=600)
            # Pass through JSON or text exactly
            try:
                out = resp.json()
                return jsonify(out), resp.status_code
            except Exception:
                return resp.text, resp.status_code, {"Content-Type": resp.headers.get("Content-Type", "text/plain")}
        except Exception as e:
            return json_error(f"Proxy to upstream failed: {e}", 502)

    # 1) If pysam is available (local machine / external server), run the real implementation
    if HAS_PYSAM:
        return _annotate_with_pysam()

    # 2) Otherwise, we are likely on Vercel serverless. Return a friendly demo result.
    return _annotate_demo()

# ----------------------- Real (pysam) path -----------------------
def _annotate_with_pysam():
    bed_file = request.files.get("file")
    gtf_bgz  = request.files.get("gtf_bgz")
    gtf_tbi  = request.files.get("gtf_tbi")
    ambiguities = request.form.get("ambiguities", "best_all")

    if not bed_file:
        return json_error("No BED file uploaded")
    if not gtf_bgz or not gtf_tbi:
        return json_error("Please upload both GTF (.gtf.bgz) and its .tbi index")

    try:
        tbx, tmpdir = _open_gtf_pair(gtf_bgz, gtf_tbi)
    except Exception as e:
        return json_error(f"GTF open failed: {e}")

    try:
        intervals = parse_bed(bed_file.read())
        rows, biotype_counter, gene_set = [], Counter(), set()

        def attr_dict(attr_str: str) -> dict:
            d = {}
            for field in attr_str.strip().split(";"):
                field = field.strip()
                if not field: continue
                if " " in field:
                    k, v = field.split(" ", 1)
                    d[k] = v.strip().strip('"')
            return d

        def overlap_len(a1,a2,b1,b2):
            return max(0, min(a2, b2) - max(a1, b1))

        def biotype_rank(bt: str) -> int:
            if not bt: return 5
            if bt == "protein_coding": return 0
            if bt.endswith("RNA"): return 2
            if bt.endswith("_decay"): return 3
            if bt.startswith("sense_"): return 4
            if bt == "antisense": return 6
            if bt.startswith("translated_"): return 7
            if bt.startswith("transcribed_"): return 8
            return 5

        def score_transcript(tx_pct, cds_pct, exon_pct, bt, has_hugo, tx_len):
            s = 3.0 * tx_pct + 2.0 * cds_pct + 1.5 * exon_pct
            s += max(0, 50 - 10 * biotype_rank(bt))
            if has_hugo: s += 5.0
            if tx_len:
                try:
                    s += min(10.0, math.log1p(tx_len)/math.log(10)*5.0)
                except: pass
            return s

        for chrom, start0, end0 in intervals:
            transcripts, exons_by_tx, cds_by_tx = {}, defaultdict(list), defaultdict(list)
            found = False

            for chrom_alt in (chrom, f"chr{chrom}"):
                try:
                    it = tbx.fetch(chrom_alt, start0, end0)
                except ValueError:
                    continue
                for line in it:
                    parts = line.rstrip("\n").split("\t")
                    if len(parts) < 9: continue
                    _seq,_src,ftype,fstart,fend,_s,strand,_ph,attrs = parts
                    fstart, fend = int(fstart), int(fend)
                    if overlap_len(start0,end0,fstart-1,fend) <= 0:
                        continue

                    a = attr_dict(attrs)
                    tx_id = a.get("transcript_id")
                    gene_id = a.get("gene_id")
                    gene_name = a.get("gene_name")
                    biotype = a.get("transcript_biotype") or a.get("gene_biotype")

                    if ftype=="transcript" and tx_id:
                        found = True
                        transcripts[tx_id] = {
                            "tx_id":tx_id, "gene_id":gene_id,"gene_name":gene_name,
                            "biotype":biotype,"strand":strand,
                            "tstart":fstart,"tend":fend
                        }
                    elif ftype=="exon" and tx_id: exons_by_tx[tx_id].append((fstart,fend))
                    elif ftype=="CDS" and tx_id:  cds_by_tx[tx_id].append((fstart,fend))
                if found: break

            candidates=[]
            for tx_id,info in transcripts.items():
                tstart,tend = info["tstart"], info["tend"]
                tx_total = max(0,tend-tstart+1)
                tx_ol = overlap_len(start0,end0,tstart-1,tend)
                tx_pct = 100*tx_ol/tx_total if tx_total>0 else 0

                ex_parts = exons_by_tx.get(tx_id,[])
                ex_cov = sum(overlap_len(start0,end0,e1-1,e2) for e1,e2 in ex_parts)
                exon_total = sum(e2-e1+1 for e1,e2 in ex_parts)
                exon_pct = 100*ex_cov/exon_total if exon_total>0 else 0

                cds_parts = cds_by_tx.get(tx_id,[])
                cds_cov = sum(overlap_len(start0,end0,c1-1,c2) for c1,c2 in cds_parts)
                cds_total = sum(c2-c1+1 for c1,c2 in cds_parts)
                cds_pct = 100*cds_cov/cds_total if cds_total>0 else 0

                score = score_transcript(tx_pct, cds_pct, exon_pct, info["biotype"], bool(info["gene_name"]), tx_total)

                candidates.append({
                    "input_chr":chrom,"input_start":start0,"input_end":end0,
                    "gene":info["gene_id"],"strand":info["strand"],"feature_biotype":info["biotype"],
                    "ensembl_id":tx_id,"hugo":info["gene_name"],
                    "tx_overlap_pct":round(tx_pct,3),"exon_overlap_pct":round(exon_pct,3),
                    "cds_overlap_pct":round(cds_pct,3),"priority_score":round(score,3)
                })

            if not candidates:
                rows.append({"input_chr":chrom,"input_start":start0,"input_end":end0,
                             "gene":None,"strand":None,"feature_biotype":None,
                             "ensembl_id":None,"hugo":None,
                             "tx_overlap_pct":0,"exon_overlap_pct":0,"cds_overlap_pct":0,
                             "priority_score":0})
                continue

            candidates.sort(key=lambda x:x["priority_score"],reverse=True)
            # “best_all” default
            top=candidates[0]["priority_score"]
            picks = candidates if request.form.get("ambiguities","best_all")=="all" else \
                    [candidates[0]] if request.form.get("ambiguities")=="best_one" else \
                    [c for c in candidates if abs(c["priority_score"]-top)<1e-6]

            for c in picks:
                rows.append(c)
                if c["feature_biotype"]: biotype_counter[c["feature_biotype"]] += 1
                if c["hugo"]: gene_set.add(c["hugo"])
                elif c["gene"]: gene_set.add(c["gene"])

        outpath = os.path.join(tempfile.gettempdir(), "annotations_upload_gtf.csv")
        pd.DataFrame(rows).to_csv(outpath, index=False)
        summary = {"unique_genes": len(gene_set), "by_biotype": dict(biotype_counter)}
        return jsonify({"rows": rows, "summary": summary, "csv_download_path": outpath})

    finally:
        try: shutil.rmtree(tmpdir, ignore_errors=True)
        except: pass

def _open_gtf_pair(gtf_file, tbi_file):
    tmpdir = tempfile.mkdtemp(prefix="gtf_")
    base = f"user_{uuid.uuid4().hex}.gtf.bgz"
    gtf_path = os.path.join(tmpdir, base)
    tbi_path = gtf_path + ".tbi"
    gtf_file.save(gtf_path)
    tbi_file.save(tbi_path)
    tbx = pysam.TabixFile(gtf_path)
    return tbx, tmpdir

# ----------------------- Demo fallback (no pysam) -----------------------
def _annotate_demo():
    bed = request.files.get("file")
    if not bed:
        return json_error("No BED file uploaded")
    intervals = parse_bed(bed.read())
    # fabricate stable fake annotations (UI demo)
    rows=[]; biotypes=["protein_coding","lincRNA","antisense","sense_intronic","processed_transcript"]
    for i,(c,s,e) in enumerate(intervals):
        rows.append({
            "input_chr": c, "input_start": s, "input_end": e,
            "gene": f"ENSG_FAKE_{i}", "strand": "+" if i%2==0 else "-",
            "feature_biotype": biotypes[i % len(biotypes)],
            "ensembl_id": f"ENST_FAKE_{i}", "hugo": f"MOCK{i%7}",
            "tx_overlap_pct": 80.0 - (i%5)*3, "exon_overlap_pct": 45.0, "cds_overlap_pct": 20.0,
            "priority_score": 99.0 - i,
        })
    summary = {
        "unique_genes": len(set(r["hugo"] for r in rows)),
        "by_biotype": dict(Counter(r["feature_biotype"] for r in rows))
    }
    outpath = os.path.join(tempfile.gettempdir(), "annotations_demo.csv")
    pd.DataFrame(rows).to_csv(outpath, index=False)
    return jsonify({"rows": rows, "summary": summary, "csv_download_path": outpath})