<p align="center">
  <h1 align="center">annota-bed</h1>
</p>

<p align="center">
  <em>Fast and easy annotation of BED files for contextual genomic analysis</em>
</p>

---

## Attribution

This project was adapted from [vladsavelyev/bed_annotation](https://github.com/vladsavelyev/bed_annotation) and extended into a full-stack tool with a Next.js frontend and Flask backend.

---

## Overview

This tool was created as a fun project to solve a recurring problem I faced:  
**Annotating genomic motifs from various reference genomes directly from a simple BED file**.  

Instead of manually intersecting BED files with GTFs and parsing outputs, `annota-bed` provides an interactive web interface where you can upload:

- A **BED file** (`chr`, `start`, `end`)
- A **reference genome GTF** (Ensembl, bgzip + tabix formatted)
- The corresponding **tabix index (.tbi)**

The tool then returns **gene- and transcript-level annotations**, interactive plots, and a downloadable CSV.

---

## How It Works

1. **User Input**
   - Upload a BED file with 3 columns: chromosome, start, end.
   - Upload a tabix-indexed **GTF** (`.gtf.bgz`) and its **.tbi index** for the chosen reference genome.
   - Ensembl-formatted GTFs for hg19 and hg38 are included in this repo for convenience.

2. **Backend (Flask)**
   - The backend uses [pysam](https://pysam.readthedocs.io/en/latest/) to query the GTF in genomic intervals.
   - Each BED region is intersected against transcripts, exons, and CDS features.

3. **Annotation Priority**
   If multiple overlapping transcripts/features are found, the tool applies a **priority system**:
   - Overlap % with transcript (higher is better)
   - Overlap % with CDS
   - Overlap % with exons
   - Biotype (protein_coding > others > \*RNA > \*_decay > sense\_* > antisense > translated\_* > transcribed\_*)
   - Presence of HUGO gene symbol (preferred over only Ensembl IDs)
   - Transcript size (longer preferred)

   The result is the "best" annotation for each region.

4. **Frontend (Next.js + Recharts)**
   - Plots distribution of annotated biotypes
   - Plots distribution of regions across chromosomes
   - Displays example gene names
   - Interactive results table with scrolling
   - Downloadable CSV with full results

---

## User Options

- **Ambiguities**  
  - `best_all` (default): Keeps all top matches if tied by priority.  
  - `best_one`: Returns only the single best match.  
  - `all`: Returns all overlapping annotations (no filtering).

- **Custom GTF Uploads**  
  Users can swap in any **Ensembl-style GTF** (properly sorted and bgzip-tabix indexed). This makes the tool flexible for other organisms and genome builds.

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
pnpm dev'''

Open http://localhost:3000 in your browser.
The Flask API will run at http://127.0.0.1:5328.