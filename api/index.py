# api/index.py
import bisect
import csv
import io
import math
import os
import random
import shutil
import tempfile
import uuid
from collections import Counter, defaultdict

from flask import Flask, jsonify, request, send_file

app = Flask(__name__)

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA_DIR = os.path.join(BASE_DIR, "data")
CCRE_PATH = os.path.join(BASE_DIR, "GRCh38-cCREs.bed")

HG38_CHROM_SIZES = {
    "1": 248956422,
    "2": 242193529,
    "3": 198295559,
    "4": 190214555,
    "5": 181538259,
    "6": 170805979,
    "7": 159345973,
    "8": 145138636,
    "9": 138394717,
    "10": 133797422,
    "11": 135086622,
    "12": 133275309,
    "13": 114364328,
    "14": 107043718,
    "15": 101991189,
    "16": 90338345,
    "17": 83257441,
    "18": 80373285,
    "19": 58617616,
    "20": 64444167,
    "21": 46709983,
    "22": 50818468,
    "X": 156040895,
    "Y": 57227415,
    "M": 16569,
}

REFERENCE_GENOMES = {
    "hg38": {
        "label": "GRCh38 / hg38",
        "description": "Ensembl-style GRCh38 gene annotations",
        "gtf": os.path.join(DATA_DIR, "hg38.sorted.gtf.bgz"),
    },
    "hg37": {
        "label": "GRCh37 / hg19",
        "description": "Ensembl GRCh37 gene annotations",
        "gtf": os.path.join(DATA_DIR, "hg37.sorted.gtf.bgz"),
    },
    "chm13": {
        "label": "T2T-CHM13v2.0",
        "description": "T2T-CHM13 gene annotations",
        "gtf": os.path.join(DATA_DIR, "chm13.sorted.gtf.bgz"),
    },
}

for ref in REFERENCE_GENOMES.values():
    ref["tbi"] = ref["gtf"] + ".tbi"

# Vercel/serverless cannot reliably ship and query these native/data-heavy files.
IS_VERCEL = os.getenv("VERCEL") == "1"

HAS_PYSAM = False
if not IS_VERCEL:
    try:
        import pysam

        HAS_PYSAM = True
    except Exception:
        HAS_PYSAM = False

_CCRE_INDEX = None

CSV_FIELDS = [
    "region_id",
    "input_chr",
    "input_start",
    "input_end",
    "region_size",
    "gene",
    "strand",
    "feature_biotype",
    "ensembl_id",
    "hugo",
    "tx_overlap_pct",
    "exon_overlap_pct",
    "cds_overlap_pct",
    "priority_score",
    "ccre_overlap_count",
    "ccre_best_type",
    "ccre_types",
    "ccre_ids",
    "ccre_max_overlap_bp",
]


# ----------------------- Shared helpers -----------------------
def normalize_chrom(chrom: str) -> str:
    chrom = (chrom or "").strip()
    if chrom.lower().startswith("chr"):
        chrom = chrom[3:]
    return chrom


def parse_bed(file_bytes: bytes):
    intervals = []
    for line in io.BytesIO(file_bytes).read().decode(errors="ignore").splitlines():
        if not line or line.startswith(("track", "browser", "#")):
            continue
        parts = line.split()
        if len(parts) < 3:
            continue
        chrom = normalize_chrom(parts[0])
        try:
            start = int(parts[1])
            end = int(parts[2])
        except ValueError:
            continue
        if end > start:
            intervals.append(
                {
                    "region_id": len(intervals) + 1,
                    "chrom": chrom,
                    "start": start,
                    "end": end,
                    "size": end - start,
                }
            )
    return intervals


def json_error(msg, code=400):
    return jsonify({"error": msg}), code


def truthy(value, default=False):
    if value is None:
        return default
    return str(value).strip().lower() in {"1", "true", "yes", "y", "on"}


def requested_ccre_mode():
    mode = request.form.get("ccre_mode")
    if mode in {"off", "unannotated", "all"}:
        return mode
    return "all" if truthy(request.form.get("include_ccre"), default=False) else "off"


def requested_ccre_permutation():
    run = truthy(request.form.get("ccre_permutation"), default=False)
    try:
        permutations = int(request.form.get("ccre_permutations", "1000"))
    except ValueError:
        permutations = 1000
    return run, max(1, min(permutations, 5000))


def safe_pct(numerator, denominator):
    return round(100.0 * numerator / denominator, 3) if denominator else 0


def mean(values):
    values = [v for v in values if isinstance(v, (int, float))]
    return round(sum(values) / len(values), 3) if values else 0


def attr_dict(attr_str: str) -> dict:
    attrs = {}
    for field in attr_str.strip().split(";"):
        field = field.strip()
        if not field:
            continue
        if " " in field:
            key, value = field.split(" ", 1)
            attrs[key] = value.strip().strip('"')
        elif "=" in field:
            key, value = field.split("=", 1)
            attrs[key] = value.strip().strip('"')
    return attrs


def overlap_len(a1, a2, b1, b2):
    return max(0, min(a2, b2) - max(a1, b1))


def biotype_rank(biotype: str) -> int:
    if not biotype:
        return 5
    if biotype == "protein_coding":
        return 0
    if biotype.endswith("RNA"):
        return 2
    if biotype.endswith("_decay"):
        return 3
    if biotype.startswith("sense_"):
        return 4
    if biotype == "antisense":
        return 6
    if biotype.startswith("translated_"):
        return 7
    if biotype.startswith("transcribed_"):
        return 8
    return 5


def score_transcript(tx_pct, cds_pct, exon_pct, biotype, has_hugo, tx_len):
    score = 3.0 * tx_pct + 2.0 * cds_pct + 1.5 * exon_pct
    score += max(0, 50 - 10 * biotype_rank(biotype))
    if has_hugo:
        score += 5.0
    if tx_len:
        try:
            score += min(10.0, math.log1p(tx_len) / math.log(10) * 5.0)
        except Exception:
            pass
    return score


def reference_payload():
    refs = []
    for key, ref in REFERENCE_GENOMES.items():
        refs.append(
            {
                "key": key,
                "label": ref["label"],
                "description": ref["description"],
                "available": os.path.exists(ref["gtf"]) and os.path.exists(ref["tbi"]),
                "gtf": os.path.basename(ref["gtf"]),
                "supports_ccre": key == "hg38" and os.path.exists(CCRE_PATH),
            }
        )
    return refs


# ----------------------- Health -----------------------
@app.get("/health")
def health():
    if IS_VERCEL:
        mode = "serverless-demo"
    else:
        mode = "full-pysam" if HAS_PYSAM else "demo-fallback-local"
    return {
        "status": "ok",
        "mode": mode,
        "references": reference_payload(),
        "ccre_available": os.path.exists(CCRE_PATH),
    }


@app.get("/references")
def references():
    return jsonify({"references": reference_payload()})


# ----------------------- Annotate -----------------------
@app.post("/annotate")
def annotate():
    """
    Local/full mode:
      - BED is uploaded by the user.
      - GTF references are read from data/ by genome selection by default.
      - Custom GTF/TBI uploads remain available for non-packaged references.

    Vercel/serverless mode:
      - Returns deterministic demo annotations so the UI still works.
    """
    if IS_VERCEL or not HAS_PYSAM:
        return _annotate_demo()
    return _annotate_with_pysam()


# ----------------------- cCRE helpers -----------------------
def _empty_ccre_summary(mode="off", enabled=False, available=False, reason=None):
    return {
        "mode": mode,
        "enabled": enabled,
        "available": available,
        "reason": reason,
        "tested_regions": 0,
        "overlap_regions": 0,
        "overlap_region_pct": 0,
        "input_overlap_region_pct": 0,
        "total_overlaps": 0,
        "by_type": {},
        "by_type_regions": {},
        "permutation": None,
        "enrichment": [],
    }


def _load_ccre_index():
    global _CCRE_INDEX
    if _CCRE_INDEX is not None:
        return _CCRE_INDEX
    if not os.path.exists(CCRE_PATH):
        return None

    by_chrom = defaultdict(list)
    catalog_counts = Counter()
    total = 0

    with open(CCRE_PATH, "r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            try:
                start = int(parts[1])
                end = int(parts[2])
            except ValueError:
                continue
            if end <= start:
                continue

            chrom = normalize_chrom(parts[0])
            ccre_id = parts[3] if len(parts) > 3 else ""
            encode_id = parts[4] if len(parts) > 4 else ""
            ccre_type = parts[5] if len(parts) > 5 else "cCRE"
            by_chrom[chrom].append((start, end, ccre_type, ccre_id, encode_id))
            catalog_counts[ccre_type] += 1
            total += 1

    indexed = {}
    for chrom, items in by_chrom.items():
        items.sort(key=lambda item: item[0])
        indexed[chrom] = {
            "items": items,
            "starts": [item[0] for item in items],
            "max_len": max((item[1] - item[0] for item in items), default=0),
        }

    _CCRE_INDEX = {
        "by_chrom": indexed,
        "catalog_counts": dict(catalog_counts),
        "catalog_total": total,
    }
    return _CCRE_INDEX


def _query_ccres(chrom, start, end):
    index = _load_ccre_index()
    if not index:
        return []

    chrom_index = index["by_chrom"].get(normalize_chrom(chrom))
    if not chrom_index:
        return []

    starts = chrom_index["starts"]
    items = chrom_index["items"]
    max_len = chrom_index["max_len"]
    lo = bisect.bisect_left(starts, max(0, start - max_len))
    hi = bisect.bisect_left(starts, end)

    overlaps = []
    for item_start, item_end, ccre_type, ccre_id, encode_id in items[lo:hi]:
        bp = overlap_len(start, end, item_start, item_end)
        if bp <= 0:
            continue
        overlaps.append(
            {
                "start": item_start,
                "end": item_end,
                "type": ccre_type,
                "ccre_id": ccre_id,
                "encode_id": encode_id,
                "overlap_bp": bp,
            }
        )
    overlaps.sort(key=lambda item: item["overlap_bp"], reverse=True)
    return overlaps


def _annotate_ccres(
    intervals,
    genome,
    ccre_mode,
    eligible_region_ids=None,
    run_permutation=False,
    num_permutations=1000,
    seed=42,
):
    if ccre_mode == "off":
        return {}, _empty_ccre_summary("off", False, os.path.exists(CCRE_PATH), "cCRE analysis disabled")
    if genome != "hg38":
        return {}, _empty_ccre_summary(ccre_mode, True, os.path.exists(CCRE_PATH), "cCRE catalog is GRCh38 only")
    if not os.path.exists(CCRE_PATH):
        return {}, _empty_ccre_summary(ccre_mode, True, False, "GRCh38 cCRE BED file not found")

    by_region = {}
    type_counts = Counter()
    region_type_counts = Counter()
    eligible_region_ids = set(eligible_region_ids) if eligible_region_ids is not None else None
    tested_intervals = [
        interval
        for interval in intervals
        if eligible_region_ids is None or interval["region_id"] in eligible_region_ids
    ]

    for interval in tested_intervals:
        overlaps = _query_ccres(interval["chrom"], interval["start"], interval["end"])
        if not overlaps:
            continue
        by_region[interval["region_id"]] = overlaps
        region_types = set()
        for overlap in overlaps:
            type_counts[overlap["type"]] += 1
            region_types.add(overlap["type"])
        for ccre_type in region_types:
            region_type_counts[ccre_type] += 1

    index = _load_ccre_index()
    catalog_counts = Counter(index["catalog_counts"]) if index else Counter()
    catalog_total = index["catalog_total"] if index else 0
    observed_total = sum(type_counts.values())
    enrichment = []
    for ccre_type, region_count in region_type_counts.most_common():
        count = type_counts.get(ccre_type, 0)
        record_fraction = count / observed_total if observed_total else 0
        region_fraction = region_count / len(tested_intervals) if tested_intervals else 0
        catalog_fraction = catalog_counts.get(ccre_type, 0) / catalog_total if catalog_total else 0
        enrichment.append(
            {
                "type": ccre_type,
                "observed_count": count,
                "region_count": region_count,
                "observed_fraction": round(record_fraction, 5),
                "region_fraction": round(region_fraction, 5),
                "catalog_fraction": round(catalog_fraction, 5),
                "fold_enrichment": round(record_fraction / catalog_fraction, 3)
                if catalog_fraction
                else None,
                "region_fold": round(region_fraction / catalog_fraction, 3)
                if catalog_fraction
                else None,
            }
        )

    permutation = None
    if run_permutation and tested_intervals:
        permutation = _run_ccre_permutation(
            tested_intervals=tested_intervals,
            observed_region_counts=region_type_counts,
            num_permutations=num_permutations,
            seed=seed,
        )

    return by_region, {
        "mode": ccre_mode,
        "enabled": True,
        "available": True,
        "reason": None,
        "tested_regions": len(tested_intervals),
        "overlap_regions": len(by_region),
        "overlap_region_pct": safe_pct(len(by_region), len(tested_intervals)),
        "input_overlap_region_pct": safe_pct(len(by_region), len(intervals)),
        "total_overlaps": observed_total,
        "by_type": dict(type_counts),
        "by_type_regions": dict(region_type_counts),
        "permutation": permutation,
        "enrichment": enrichment,
    }


def _random_interval_for_width(width, rng, choice_cache):
    if width in choice_cache:
        cumulative, choices, total_weight = choice_cache[width]
    else:
        cumulative = []
        choices = []
        total_weight = 0
        for chrom, chrom_len in HG38_CHROM_SIZES.items():
            max_start = chrom_len - width
            if max_start < 0:
                continue
            weight = max_start + 1
            total_weight += weight
            cumulative.append(total_weight)
            choices.append((chrom, max_start))
        choice_cache[width] = (cumulative, choices, total_weight)

    if not choices or total_weight <= 0:
        return None

    draw = rng.randrange(total_weight)
    idx = bisect.bisect_right(cumulative, draw)
    chrom, max_start = choices[idx]
    start = rng.randint(0, max_start)
    return {"chrom": chrom, "start": start, "end": start + width, "size": width}

def _run_ccre_permutation(tested_intervals, observed_region_counts, num_permutations, seed):
    index = _load_ccre_index()
    if not index:
        return None

    classes = sorted(index["catalog_counts"].keys())
    observed = {ccre_class: observed_region_counts.get(ccre_class, 0) for ccre_class in classes}
    expected_sums = Counter()
    expected_gte_observed = Counter()
    rng = random.Random(seed)
    widths = [max(1, interval["end"] - interval["start"]) for interval in tested_intervals]
    choice_cache = {}

    for _ in range(num_permutations):
        perm_counts = Counter()
        for width in widths:
            random_interval = _random_interval_for_width(width, rng, choice_cache)
            if not random_interval:
                continue
            overlaps = _query_ccres(random_interval["chrom"], random_interval["start"], random_interval["end"])
            for ccre_type in {overlap["type"] for overlap in overlaps if overlap.get("type")}:
                perm_counts[ccre_type] += 1

        for ccre_class in classes:
            perm_count = perm_counts.get(ccre_class, 0)
            expected_sums[ccre_class] += perm_count
            if perm_count >= observed[ccre_class]:
                expected_gte_observed[ccre_class] += 1

    results = []
    for ccre_class in classes:
        expected = expected_sums[ccre_class] / num_permutations
        observed_count = observed[ccre_class]
        results.append(
            {
                "type": ccre_class,
                "observed": observed_count,
                "expected": round(expected, 4),
                "fold_enrichment": round(observed_count / (expected + 1e-5), 4),
                "empirical_p": round(expected_gte_observed[ccre_class] / num_permutations, 5),
            }
        )

    results.sort(key=lambda row: (row["fold_enrichment"], row["observed"]), reverse=True)
    return {
        "enabled": True,
        "method": "width-preserving genome randomization",
        "permutations": num_permutations,
        "seed": seed,
        "tested_regions": len(tested_intervals),
        "results": results,
    }


def _row_ccre_fields(region_id, ccre_by_region):
    overlaps = ccre_by_region.get(region_id, [])
    if not overlaps:
        return {
            "ccre_overlap_count": 0,
            "ccre_best_type": None,
            "ccre_types": None,
            "ccre_ids": None,
            "ccre_max_overlap_bp": 0,
        }

    types = sorted({item["type"] for item in overlaps if item["type"]})
    ids = [item["ccre_id"] or item["encode_id"] for item in overlaps[:5]]
    ids = [item for item in ids if item]
    return {
        "ccre_overlap_count": len(overlaps),
        "ccre_best_type": overlaps[0]["type"],
        "ccre_types": ", ".join(types) if types else None,
        "ccre_ids": ", ".join(ids) if ids else None,
        "ccre_max_overlap_bp": overlaps[0]["overlap_bp"],
    }


# ----------------------- Summary / output helpers -----------------------
def _build_summary(intervals, rows, ccre_by_region):
    region_count = len(intervals)
    by_chrom = Counter(interval["chrom"] for interval in intervals)
    lengths = [interval["size"] for interval in intervals]
    matched_region_ids = {
        row["region_id"]
        for row in rows
        if row.get("gene") or row.get("hugo") or row.get("ensembl_id")
    }
    genes = {
        row.get("hugo") or row.get("gene")
        for row in rows
        if row.get("hugo") or row.get("gene")
    }
    by_biotype = Counter(row.get("feature_biotype") for row in rows if row.get("feature_biotype"))
    biotype_region_counts = Counter()

    rows_by_region = defaultdict(list)
    for row in rows:
        rows_by_region[row["region_id"]].append(row)

    for region_rows in rows_by_region.values():
        for biotype in {row.get("feature_biotype") for row in region_rows if row.get("feature_biotype")}:
            biotype_region_counts[biotype] += 1

    region_classes = Counter()
    for interval in intervals:
        region_id = interval["region_id"]
        region_rows = rows_by_region.get(region_id, [])
        if any(float(row.get("cds_overlap_pct") or 0) > 0 for row in region_rows):
            region_classes["CDS overlap"] += 1
        elif any(float(row.get("exon_overlap_pct") or 0) > 0 for row in region_rows):
            region_classes["exonic non-CDS"] += 1
        elif any(float(row.get("tx_overlap_pct") or 0) > 0 for row in region_rows):
            region_classes["transcript/intronic"] += 1
        elif ccre_by_region.get(region_id):
            region_classes["cCRE only"] += 1
        else:
            region_classes["unannotated"] += 1

    gene_rows = [row for row in rows if row.get("gene") or row.get("hugo")]

    return {
        "total_regions": region_count,
        "reported_rows": len(rows),
        "annotated_regions": len(matched_region_ids),
        "unannotated_regions": max(0, region_count - len(matched_region_ids)),
        "annotated_region_pct": safe_pct(len(matched_region_ids), region_count),
        "unique_genes": len(genes),
        "median_region_bp": sorted(lengths)[len(lengths) // 2] if lengths else 0,
        "mean_region_bp": mean(lengths),
        "mean_tx_overlap_pct": mean([row.get("tx_overlap_pct", 0) for row in gene_rows]),
        "mean_exon_overlap_pct": mean([row.get("exon_overlap_pct", 0) for row in gene_rows]),
        "mean_cds_overlap_pct": mean([row.get("cds_overlap_pct", 0) for row in gene_rows]),
        "by_biotype": dict(by_biotype),
        "by_biotype_regions": dict(biotype_region_counts),
        "by_chrom": dict(by_chrom),
        "region_classes": dict(region_classes),
    }


def _write_csv(rows):
    outpath = os.path.join(tempfile.gettempdir(), f"annotations_{uuid.uuid4().hex}.csv")
    with open(outpath, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=CSV_FIELDS)
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key) for key in CSV_FIELDS})
    return outpath


def _response_payload(rows, intervals, ccre_by_region, ccre_summary, csv_path, genome, source):
    ref = REFERENCE_GENOMES.get(genome, {})
    return jsonify(
        {
            "rows": rows,
            "summary": _build_summary(intervals, rows, ccre_by_region),
            "ccre": ccre_summary,
            "csv_download_path": csv_path,
            "reference": {
                "genome": genome,
                "label": ref.get("label", "Custom reference"),
                "source": source,
            },
        }
    )


# ----------------------- DEMO (serverless-safe) -----------------------
def _annotate_demo():
    bed = request.files.get("file")
    if not bed:
        return json_error("No BED file uploaded")

    genome = request.form.get("genome", "hg38")
    ccre_mode = requested_ccre_mode()
    run_ccre_permutation, ccre_permutations = requested_ccre_permutation()
    intervals = parse_bed(bed.read())
    ccre_by_region, ccre_summary = ({}, _empty_ccre_summary(ccre_mode, ccre_mode != "off", False, "demo mode"))

    rows = []
    biotypes = ["protein_coding", "lincRNA", "antisense", "sense_intronic", "processed_transcript"]
    for idx, interval in enumerate(intervals):
        rows.append(
            {
                "region_id": interval["region_id"],
                "input_chr": interval["chrom"],
                "input_start": interval["start"],
                "input_end": interval["end"],
                "region_size": interval["size"],
                "gene": f"ENSG_FAKE_{idx}",
                "strand": "+" if idx % 2 == 0 else "-",
                "feature_biotype": biotypes[idx % len(biotypes)],
                "ensembl_id": f"ENST_FAKE_{idx}",
                "hugo": f"MOCK{idx % 7}",
                "tx_overlap_pct": round(80.0 - (idx % 5) * 3, 3),
                "exon_overlap_pct": 45.0,
                "cds_overlap_pct": 20.0,
                "priority_score": round(99.0 - idx, 3),
            }
        )

    if not IS_VERCEL:
        eligible_region_ids = None
        if ccre_mode == "unannotated":
            eligible_region_ids = []
        ccre_by_region, ccre_summary = _annotate_ccres(
            intervals,
            genome,
            ccre_mode,
            eligible_region_ids,
            run_ccre_permutation,
            ccre_permutations,
        )
    for row in rows:
        row.update(_row_ccre_fields(row["region_id"], ccre_by_region))

    outpath = _write_csv(rows)
    return _response_payload(rows, intervals, ccre_by_region, ccre_summary, outpath, genome, "demo")


# ----------------------- REAL (local pysam) -----------------------
def _annotate_with_pysam():
    bed_file = request.files.get("file")
    if not bed_file:
        return json_error("No BED file uploaded")

    genome = request.form.get("genome", "hg38")
    ambiguities = request.form.get("ambiguities", "best_all")
    ccre_mode = requested_ccre_mode()
    run_ccre_permutation, ccre_permutations = requested_ccre_permutation()
    reference_mode = request.form.get("reference_mode", "builtin")

    tbx = None
    tmpdir = None
    source = "built-in"
    try:
        if reference_mode == "custom" or request.files.get("gtf_bgz"):
            gtf_bgz = request.files.get("gtf_bgz")
            gtf_tbi = request.files.get("gtf_tbi")
            if not gtf_bgz or not gtf_tbi:
                return json_error("Please upload both GTF (.gtf.bgz) and its .tbi index")
            tbx, tmpdir = _open_gtf_pair(gtf_bgz, gtf_tbi)
            source = "custom upload"
        else:
            tbx = _open_builtin_reference(genome)

        intervals = parse_bed(bed_file.read())
        rows = []

        for interval in intervals:
            chrom = interval["chrom"]
            start0 = interval["start"]
            end0 = interval["end"]
            transcripts = {}
            exons_by_tx = defaultdict(list)
            cds_by_tx = defaultdict(list)
            found_transcript = False

            for chrom_alt in (chrom, f"chr{chrom}"):
                try:
                    iterator = tbx.fetch(chrom_alt, start0, end0)
                except ValueError:
                    continue
                for line in iterator:
                    parts = line.rstrip("\n").split("\t")
                    if len(parts) < 9:
                        continue
                    _seq, _src, feature_type, fstart, fend, _score, strand, _phase, attrs = parts
                    try:
                        fstart = int(fstart)
                        fend = int(fend)
                    except ValueError:
                        continue
                    if overlap_len(start0, end0, fstart - 1, fend) <= 0:
                        continue

                    parsed_attrs = attr_dict(attrs)
                    tx_id = parsed_attrs.get("transcript_id") or parsed_attrs.get("transcript")
                    gene_id = parsed_attrs.get("gene_id") or parsed_attrs.get("gene")
                    gene_name = parsed_attrs.get("gene_name") or parsed_attrs.get("Name")
                    biotype = (
                        parsed_attrs.get("transcript_biotype")
                        or parsed_attrs.get("transcript_type")
                        or parsed_attrs.get("gene_biotype")
                        or parsed_attrs.get("gene_type")
                    )

                    if feature_type == "transcript" and tx_id:
                        found_transcript = True
                        transcripts[tx_id] = {
                            "tx_id": tx_id,
                            "gene_id": gene_id,
                            "gene_name": gene_name,
                            "biotype": biotype,
                            "strand": strand,
                            "tstart": fstart,
                            "tend": fend,
                        }
                    elif feature_type == "exon" and tx_id:
                        exons_by_tx[tx_id].append((fstart, fend))
                    elif feature_type == "CDS" and tx_id:
                        cds_by_tx[tx_id].append((fstart, fend))
                if found_transcript:
                    break

            candidates = []
            for tx_id, info in transcripts.items():
                tstart, tend = info["tstart"], info["tend"]
                tx_total = max(0, tend - tstart + 1)
                tx_ol = overlap_len(start0, end0, tstart - 1, tend)
                tx_pct = safe_pct(tx_ol, tx_total)

                exons = exons_by_tx.get(tx_id, [])
                exon_cov = sum(overlap_len(start0, end0, exon_start - 1, exon_end) for exon_start, exon_end in exons)
                exon_total = sum(exon_end - exon_start + 1 for exon_start, exon_end in exons)
                exon_pct = safe_pct(exon_cov, exon_total)

                cds_parts = cds_by_tx.get(tx_id, [])
                cds_cov = sum(overlap_len(start0, end0, cds_start - 1, cds_end) for cds_start, cds_end in cds_parts)
                cds_total = sum(cds_end - cds_start + 1 for cds_start, cds_end in cds_parts)
                cds_pct = safe_pct(cds_cov, cds_total)

                score = score_transcript(
                    tx_pct,
                    cds_pct,
                    exon_pct,
                    info["biotype"],
                    bool(info["gene_name"]),
                    tx_total,
                )

                candidates.append(
                    {
                        "region_id": interval["region_id"],
                        "input_chr": chrom,
                        "input_start": start0,
                        "input_end": end0,
                        "region_size": interval["size"],
                        "gene": info["gene_id"],
                        "strand": info["strand"],
                        "feature_biotype": info["biotype"],
                        "ensembl_id": tx_id,
                        "hugo": info["gene_name"],
                        "tx_overlap_pct": tx_pct,
                        "exon_overlap_pct": exon_pct,
                        "cds_overlap_pct": cds_pct,
                        "priority_score": round(score, 3),
                    }
                )

            if not candidates:
                rows.append(
                    {
                        "region_id": interval["region_id"],
                        "input_chr": chrom,
                        "input_start": start0,
                        "input_end": end0,
                        "region_size": interval["size"],
                        "gene": None,
                        "strand": None,
                        "feature_biotype": None,
                        "ensembl_id": None,
                        "hugo": None,
                        "tx_overlap_pct": 0,
                        "exon_overlap_pct": 0,
                        "cds_overlap_pct": 0,
                        "priority_score": 0,
                    }
                )
                continue

            candidates.sort(key=lambda row: row["priority_score"], reverse=True)
            if ambiguities == "all":
                picks = candidates
            elif ambiguities == "best_one":
                picks = [candidates[0]]
            else:
                top = candidates[0]["priority_score"]
                picks = [row for row in candidates if abs(row["priority_score"] - top) < 1e-6]
            rows.extend(picks)

        matched_region_ids = {
            row["region_id"]
            for row in rows
            if row.get("gene") or row.get("hugo") or row.get("ensembl_id")
        }
        eligible_region_ids = None
        if ccre_mode == "unannotated":
            eligible_region_ids = [
                interval["region_id"]
                for interval in intervals
                if interval["region_id"] not in matched_region_ids
            ]
        ccre_by_region, ccre_summary = _annotate_ccres(
            intervals,
            genome,
            ccre_mode,
            eligible_region_ids,
            run_ccre_permutation,
            ccre_permutations,
        )
        for row in rows:
            row.update(_row_ccre_fields(row["region_id"], ccre_by_region))

        outpath = _write_csv(rows)
        return _response_payload(rows, intervals, ccre_by_region, ccre_summary, outpath, genome, source)

    except Exception as exc:
        return json_error(f"Annotation failed: {exc}", 500)
    finally:
        try:
            if tbx:
                tbx.close()
        except Exception:
            pass
        if tmpdir:
            shutil.rmtree(tmpdir, ignore_errors=True)


def _open_builtin_reference(genome):
    if genome not in REFERENCE_GENOMES:
        raise ValueError(f"Unknown reference genome: {genome}")
    ref = REFERENCE_GENOMES[genome]
    if not os.path.exists(ref["gtf"]) or not os.path.exists(ref["tbi"]):
        raise FileNotFoundError(f"Missing built-in reference files for {ref['label']}")
    return pysam.TabixFile(ref["gtf"])


def _open_gtf_pair(gtf_file, tbi_file):
    tmpdir = tempfile.mkdtemp(prefix="gtf_")
    base = f"user_{uuid.uuid4().hex}.gtf.bgz"
    gtf_path = os.path.join(tmpdir, base)
    tbi_path = gtf_path + ".tbi"
    gtf_file.save(gtf_path)
    tbi_file.save(tbi_path)
    tbx = pysam.TabixFile(gtf_path)
    return tbx, tmpdir


# ----------------------- Download -----------------------
@app.get("/download")
def download():
    path = request.args.get("path")
    if not path:
        return jsonify({"error": "file not found"}), 404

    abs_path = os.path.abspath(path)
    tmp_root = os.path.abspath(tempfile.gettempdir())
    if not abs_path.startswith(tmp_root + os.sep):
        return jsonify({"error": "invalid download path"}), 400
    if not os.path.exists(abs_path):
        return jsonify({"error": "file not found"}), 404
    return send_file(abs_path, as_attachment=True, download_name="annotations.csv")
