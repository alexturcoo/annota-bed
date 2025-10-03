<p align="center">
  <h1 align="center">annota-bed</h1>
</p>

<p align="center">
  <em>Fast and easy annotation of BED files for contextual genomic analysis</em>
</p>

<p align="center">
  <img src="imgs/img1.png" alt="annota-bed demo" width="700">
</p>

---

## Attribution

This project was inspired by [vladsavelyev/bed_annotation](https://github.com/vladsavelyev/bed_annotation) and extended into a full-stack tool with a Next.js frontend and Flask backend.

## Overview

`annota-bed` is a lightweight full-stack tool for annotating BED files against various reference genome annotations.  

Instead of manually intersecting BEDs with GTFs and parsing outputs, `annota-bed` provides an interactive web interface where you can upload:

- A **BED file** (`chr`, `start`, `end`)  
- A **reference genome GTF** (gtf format, bgzip + tabix indexed)  
- The corresponding **tabix index (.tbi)**  

The tool then returns **gene- and transcript-level annotations**, interactive plots, and a downloadable CSV.  

---

## Supported Reference Genomes

We currently support three reference genomes out-of-the-box, all bgzipped and tabix-indexed for fast queries:  

| Genome Build | Source | Notes |
|--------------|--------|-------|
| **GRCh37 / hg19** | [Ensembl GRCh37 release-115 GTF](https://ftp.ensembl.org/pub/grch37/release-115/gtf/homo_sapiens/) | Standard Ensembl release |
| **GRCh38 / hg38** | [Ensembl current GTF](https://ftp.ensembl.org/pub/current_gtf/homo_sapiens/) | Standard Ensembl release |
| **T2T-CHM13v2.0** | [HPRC CHM13 GFF3 annotation](https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/annotation/chm13.draft_v2.0.gene_annotation.gff3) | Converted from GFF3 → GTF, then bgzipped and tabix-indexed |

Users can seamlessly switch between hg19, hg38, and CHM13 when annotating motifs, enabling cross-genome comparisons.  

---

## How It Works

1. **User Input**
   - Upload a BED file with 3 columns: chromosome, start, end.  
   - Upload a tabix-indexed GTF (`.gtf.bgz`) and its `.tbi` index for the chosen reference genome.  
   - Pre-packaged GTFs for hg19, hg38, and CHM13 are included.  

2. **Backend (Flask)**
   - Uses [pysam](https://pysam.readthedocs.io/en/latest/) to query GTF annotations in genomic intervals.  
   - Each BED region is intersected against transcripts, exons, and CDS features.  

3. **Annotation Priority**
   If multiple overlapping transcripts/features are found, the tool applies a priority system:
   - Transcript overlap %  
   - CDS overlap %  
   - Exon overlap %  
   - Biotype (protein_coding > others > \*RNA > \*_decay > sense\_* > antisense > translated\_* > transcribed\_*)  
   - Presence of HUGO gene symbol  
   - Transcript length  

4. **Frontend (Next.js + Recharts)**
   - Plots annotated biotypes  
   - Chromosome-wise distribution plots  
   - Displays gene names and IDs  
   - Interactive results table  
   - Exportable CSV  

---

## User Options

- **Ambiguities**
  - `best_all` (default): Keeps all top matches if tied by priority.  
  - `best_one`: Returns only the single best match.  
  - `all`: Returns all overlapping annotations (no filtering).  

- **Custom GTF Uploads**  
  Users can upload any **Ensembl-style GTF** (properly sorted, bgzipped, tabix-indexed) for other organisms/genomes.  

---

## Running Locally

Clone the repo:

```bash
git clone https://github.com/<your-username>/annota-bed.git
cd annota-bed

# Frontend
pnpm install
# or npm install / yarn install

# Backend (Python)
pip install -r requirements.txt

# Start Flask backend
pnpm run flask-dev

# Start Next.js frontend
pnpm dev