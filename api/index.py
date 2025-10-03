# api/index.py
from flask import Flask, request, jsonify, send_file
import os, io, math, tempfile
from collections import defaultdict, Counter
import pandas as pd
import pysam  # ensure pysam is in requirements.txt

app = Flask(__name__)

# ---------- CONFIG: point at /data ----------
BASE_DIR = os.path.dirname(os.path.dirname(__file__))  # project root
DATA_DIR = os.path.join(BASE_DIR, "data")

GTF_PATHS = {
    "hg38": os.path.join(DATA_DIR, "hg38.sorted.gtf.bgz"),  # GRCh38
    "hg19": os.path.join(DATA_DIR, "hg37.sorted.gtf.bgz"),  # GRCh37/hg19
    # Add CHM13 when ready:
    # "chm13": os.path.join(DATA_DIR, "chm13.sorted.gtf.bgz"),
}

# ---------- helpers ----------
def parse_bed(file_bytes: bytes):
    """Return list of (chrom, start0, end0) from BED (0-based, half-open)."""
    out = []
    for line in io.BytesIO(file_bytes).read().decode(errors="ignore").splitlines():
        line = line.strip()
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

def attr_dict(attr_str: str) -> dict:
    """Parse GTF attributes column to dict."""
    d = {}
    for field in attr_str.strip().split(";"):
        field = field.strip()
        if not field:
            continue
        if " " in field:
            k, v = field.split(" ", 1)
            d[k] = v.strip().strip('"')
    return d

def overlap_len(a1,a2,b1,b2):
    """Overlap length for half-open intervals [a1,a2) and [b1,b2)."""
    left = max(a1, b1)
    right = min(a2, b2)
    return max(0, right - left)

# ---- priority components (your rules) ----
def biotype_rank(bt: str) -> int:
    # protein_coding > others > *RNA > *_decay > sense_* > antisense > translated_* > transcribed_*
    if not bt: return 5
    if bt == "protein_coding": return 0
    if bt.endswith("RNA"): return 2
    if bt.endswith("_decay"): return 3
    if bt.startswith("sense_"): return 4
    if bt == "antisense": return 6
    if bt.startswith("translated_"): return 7
    if bt.startswith("transcribed_"): return 8
    return 5

def tsl_rank(tsl):
    # 1 > NA > others > 2 > 3 > 4 > 5
    if tsl in ("1", 1): return 0
    if tsl in (None, "NA", "Not_available"): return 1
    if tsl in ("2", 2): return 3
    if tsl in ("3", 3): return 4
    if tsl in ("4", 4): return 5
    if tsl in ("5", 5): return 6
    return 2

def score_transcript(tx_pct, cds_pct, exon_pct, bt, tsl, has_hugo, tx_len):
    # Heavier weight on transcript %, then CDS, then exons
    s = 0.0
    s += 3.0 * tx_pct
    s += 2.0 * cds_pct
    s += 1.5 * exon_pct
    s += max(0, 50 - 10 * biotype_rank(bt))
    s += max(0, 20 - 3 * tsl_rank(tsl))
    if has_hugo: s += 5.0
    if tx_len:
        try:
            s += min(10.0, math.log1p(tx_len) / math.log(10) * 5.0)
        except Exception:
            pass
    return s

def open_gtf(reference: str) -> pysam.TabixFile:
    path = GTF_PATHS.get(reference)
    if not path or not os.path.exists(path):
        raise FileNotFoundError(f"GTF not found for {reference}: {path}")
    tbi = path + ".tbi"
    if not os.path.exists(tbi):
        raise FileNotFoundError(f"Tabix index (.tbi) not found for {path}")
    return pysam.TabixFile(path)

# ---------- routes ----------
@app.get("/health")
def health():
    return {"status": "ok", "data_dir": DATA_DIR, "available_refs": list(GTF_PATHS.keys())}

@app.post("/annotate")
def annotate():
    if "file" not in request.files:
        return jsonify({
            "rows": [],
            "summary": {"unique_genes": 0, "by_biotype": {}},
            "csv_download_path": None,
            "error": "No file uploaded",
        }), 400

    reference = request.form.get("reference", "hg38")
    ambiguities = request.form.get("ambiguities", "best_all")  # "best_one" | "best_all" | "all"

    try:
        tbx = open_gtf(reference)
    except Exception as e:
        return jsonify({"error": str(e)}), 500

    intervals = parse_bed(request.files["file"].read())
    if not intervals:
        return jsonify({"rows": [], "summary": {"unique_genes": 0, "by_biotype": {}}, "csv_download_path": None})

    rows = []
    biotype_counter = Counter()
    gene_set = set()

    for chrom, start0, end0 in intervals:
        # collect features for candidates in this interval
        transcripts = {}
        exons_by_tx = defaultdict(list)
        cds_by_tx   = defaultdict(list)

        # Try both chrom naming styles
        for chrom_alt in (chrom, f"chr{chrom}"):
            try:
                it = tbx.fetch(chrom_alt, start0, end0)
            except ValueError:
                continue
            for line in it:
                # GTF columns
                # seqname, source, feature, start, end, score, strand, frame, attributes
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 9:
                    continue
                _seq, _src, ftype, fstart, fend, _score, strand, _phase, attrs = parts
                fstart = int(fstart)
                fend = int(fend)

                # Convert GTF 1-based inclusive to half-open [fstart-1, fend)
                if overlap_len(start0, end0, fstart - 1, fend) <= 0:
                    continue

                a = attr_dict(attrs)
                tx_id = a.get("transcript_id")
                gene_id = a.get("gene_id")
                gene_name = a.get("gene_name") or a.get("gene")
                biotype = a.get("transcript_biotype") or a.get("gene_biotype")
                tsl = a.get("transcript_support_level") or a.get("tsl")

                if ftype == "transcript" and tx_id:
                    transcripts[tx_id] = {
                        "tx_id": tx_id,
                        "gene_id": gene_id,
                        "gene_name": gene_name,
                        "biotype": biotype,
                        "strand": strand,
                        "tstart": fstart,
                        "tend": fend,
                        "tsl": tsl,
                    }
                elif ftype == "exon" and tx_id:
                    exons_by_tx[tx_id].append((fstart, fend))
                elif ftype == "CDS" and tx_id:
                    cds_by_tx[tx_id].append((fstart, fend))
            if transcripts:
                break  # found hits with one naming scheme

        # Compute per-transcript metrics
        candidates = []
        for tx_id, info in transcripts.items():
            tstart, tend = info["tstart"], info["tend"]
            tx_total = max(0, tend - tstart + 1)
            tx_ol    = overlap_len(start0, end0, tstart - 1, tend)
            tx_pct   = 100.0 * tx_ol / tx_total if tx_total > 0 else 0.0

            ex_parts = exons_by_tx.get(tx_id, [])
            if ex_parts:
                ex_total = sum(e2 - e1 + 1 for e1, e2 in ex_parts)
                ex_cov   = sum(overlap_len(start0, end0, e1 - 1, e2) for e1, e2 in ex_parts)
                exon_pct = 100.0 * ex_cov / ex_total if ex_total > 0 else 0.0
            else:
                exon_pct = 0.0

            cds_parts = cds_by_tx.get(tx_id, [])
            if cds_parts:
                cds_total = sum(c2 - c1 + 1 for c1, c2 in cds_parts)
                cds_cov   = sum(overlap_len(start0, end0, c1 - 1, c2) for c1, c2 in cds_parts)
                cds_pct   = 100.0 * cds_cov / cds_total if cds_total > 0 else 0.0
            else:
                cds_pct = 0.0

            score = score_transcript(
                tx_pct, cds_pct, exon_pct,
                info["biotype"], info["tsl"],
                bool(info["gene_name"]), tx_total
            )

            candidates.append({
                "input_chr": chrom,
                "input_start": start0,
                "input_end": end0,
                "gene": info["gene_id"],
                "exon": None,
                "strand": info["strand"],
                "feature_biotype": info["biotype"],
                "ensembl_id": tx_id,
                "tsl": str(info["tsl"] or "NA"),
                "hugo": info["gene_name"],
                "tx_overlap_pct": round(tx_pct, 3),
                "exon_overlap_pct": round(exon_pct, 3),
                "cds_overlap_pct": round(cds_pct, 3),
                "priority_score": round(score, 3),
                "ori_fields": {"transcript_start": tstart, "transcript_end": tend},
            })

        if not candidates:
            rows.append({
                "input_chr": chrom, "input_start": start0, "input_end": end0,
                "gene": None, "exon": None, "strand": None, "feature_biotype": None,
                "ensembl_id": None, "tsl": None, "hugo": None,
                "tx_overlap_pct": 0.0, "exon_overlap_pct": 0.0, "cds_overlap_pct": 0.0,
                "priority_score": 0.0, "ori_fields": {}
            })
            continue

        candidates.sort(key=lambda x: x["priority_score"], reverse=True)

        # ambiguity resolution
        if ambiguities in ("all", "all_ask"):
            picks = candidates
        elif ambiguities in ("best_all", "best_ask"):
            top = candidates[0]["priority_score"]
            picks = [c for c in candidates if abs(c["priority_score"] - top) < 1e-6]
        else:  # best_one
            picks = [candidates[0]]

        for c in picks:
            rows.append(c)
            if c["feature_biotype"]:
                biotype_counter[c["feature_biotype"]] += 1
            if c["hugo"]:
                gene_set.add(c["hugo"])
            elif c["gene"]:
                gene_set.add(c["gene"])

    summary = {"unique_genes": len(gene_set), "by_biotype": dict(biotype_counter)}

    outpath = os.path.join(tempfile.gettempdir(), "annotations_local_gtf.csv")
    pd.DataFrame(rows).to_csv(outpath, index=False)

    return jsonify({"rows": rows, "summary": summary, "csv_download_path": outpath})

@app.get("/download")
def download():
    path = request.args.get("path")
    if not path or not os.path.exists(path):
        return jsonify({"error":"file not found"}), 404
    return send_file(path, as_attachment=True, download_name="annotations.csv")