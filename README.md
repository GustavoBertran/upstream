# Welcome to upstream

`upstream` is a command-line preprocessing pipeline for high-throughput sequencing (HTS) data.
It wraps FastQC, fastp, STAR, Salmon, Bowtie2, MACS2, and Bismark into four simple commands
and produces output matrices compatible with the [OBAMA pipeline](https://github.com/AOG-Lab/OBAMA).

---

## Installation

### Requirements

- Linux or macOS
- [Conda](https://docs.conda.io/en/latest/miniconda.html) (Miniconda or Anaconda)

### Steps

```bash
# 1. Clone the repository
git clone https://github.com/GustavoBertran/upstream.git
cd upstream

# 2. Create the conda environment (installs all bioinformatics tools)
conda env create -f environment.yml
conda activate upstream

# 3. Install upstream itself
pip install -e .

# 4. Confirm it works
upstream --help
```

---

## Web Interface

`upstream` includes a browser-based UI that lets you run pipelines without typing commands.
Launch it with:

```bash
upstream serve
```

Then open **http://localhost:8421** in your browser.

The interface lets you:
- **Browse** your filesystem to fill in file/folder paths (click the `...` button next to any input)
- **Select a track** in the sidebar (RNA-seq, ATAC-seq, Methylation, QC, Download Data)
- **Stream output** live as the pipeline runs
- **Download GEO data** — type a GSE accession, assign disease/control groups, and download FASTQ files with one click
- **Process 450K/EPIC methylation arrays** — select "Methylation → 450K / EPIC array (GEO beta matrix)" and provide a beta-value CSV + metadata CSV directly

To use a different port:

```bash
upstream serve --port 9000
```

---

## Quickstart

### 1. Create a samplesheet

Create a CSV file describing your samples. The `group` column must be either `disease`
or `control` — these labels are required by the OBAMA pipeline.

```csv
name,group,r1,r2
GSM123456,disease,/data/tumor_R1.fastq.gz,/data/tumor_R2.fastq.gz
GSM123457,control,/data/normal_R1.fastq.gz,/data/normal_R2.fastq.gz
```

### 2. Run the pipeline for your data type

**Quality control (run first for any data type):**
```bash
upstream qc \
  --samples samples.csv \
  --outdir results/qc/
```

**RNA-seq (Salmon, alignment-free — default):**
```bash
upstream rnaseq \
  --samples samples.csv \
  --aligner salmon \
  --salmon-index /path/to/salmon_hg38 \
  --outdir results/rnaseq/
```

**RNA-seq (STAR, genome alignment + gene counts):**
```bash
upstream rnaseq \
  --samples samples.csv \
  --aligner star \
  --star-index /path/to/star_hg38 \
  --outdir results/rnaseq/
```

**Output format (OBAMA vs DESeq2/edgeR/limma):**

RNA-seq emits the OBAMA matrix by default. Add `--format matrix` (or `--format both`)
to also write `counts_matrix.csv` (features × samples) + `coldata.csv` — the raw-count
inputs DESeq2, edgeR, and limma-voom expect. For gene-level Salmon counts, pass a
`--gtf` (or `--tx2gene`) so transcripts are aggregated to genes; STAR is already
gene-level. Both layouts use raw counts (Salmon `NumReads`, STAR `ReadsPerGene`).

```bash
upstream rnaseq --samples samples.csv --aligner salmon \
  --salmon-index /path/to/salmon_hg38 \
  --gtf /path/to/gencode.annotation.gtf.gz \
  --format both --outdir results/rnaseq/
```

**ATAC-seq:**
```bash
upstream atacseq \
  --samples samples.csv \
  --bowtie2-index /path/to/bowtie2/hg38 \
  --outdir results/atacseq/
```

With `--format matrix` (or `both`), ATAC-seq builds a consensus peak set across
samples and counts reads per peak (deeptools `multiBamSummary`), writing a
`counts_matrix.csv` (peaks × samples) + `coldata.csv` for **DESeq2/edgeR**
differential accessibility — instead of the default MACS2-score OBAMA matrix
(the score isn't a count). See `content/atacseq_export_formats.md`.

**ChIP-seq:**
```bash
upstream chipseq \
  --samples samples.csv \
  --bowtie2-index /path/to/bowtie2/hg38 \
  --peak-type narrow \
  --outdir results/chipseq/
```

ChIP-seq mirrors ATAC-seq but with two ChIP essentials:

- **Input control** — the samplesheet may add an optional `control` column naming
  each ChIP sample's input (by `name`); rows with `group=input` are aligned to
  provide the MACS2 `-c` background but are not peak-called or placed in the matrix.
  Input is optional (MACS2 runs without it).
- **`--peak-type`** — `narrow` (default; TFs, H3K4me3, H3K27ac) or `broad`
  (H3K27me3, H3K9me3, H3K36me3 → MACS2 `--broad`).

Default output is the **DESeq2/edgeR** consensus-peak count matrix (the validated
downstream for peak data). `--format obama` is available but **experimental** for
ChIP: OBAMA's gene-based interpretation (GO/Enrichr/STRING) needs peaks annotated
to genes first. See `content/chipseq_export_formats.md`.

Example ChIP samplesheet:
```csv
name,group,r1,r2,control
H3K27ac_tumor,disease,t_R1.fq.gz,t_R2.fq.gz,input_tumor
H3K27ac_normal,control,n_R1.fq.gz,n_R2.fq.gz,input_normal
input_tumor,input,it_R1.fq.gz,it_R2.fq.gz,
input_normal,input,in_R1.fq.gz,in_R2.fq.gz,
```

**DNA Methylation — WGBS (Bismark pipeline):**
```bash
upstream methylation \
  --method wgbs \
  --samples samples.csv \
  --bismark-genome /path/to/bismark_hg38 \
  --outdir results/methylation/
```

Like RNA-seq, methylation accepts `--format matrix` (or `both`) to also write an
`mvalues_matrix.csv` (features × samples) + `coldata.csv` for **limma** — M-values
(`log2(beta/(1-beta))`), which limma models instead of beta/percent. Works for both
`--method wgbs` and `--method array`.

**DNA Methylation — Illumina array (450K / EPIC, from GEO):**

No FASTQ files needed — provide the beta-value matrix downloaded from GEO
and a metadata CSV with `geo_accession` and `disease.state` columns:

```bash
upstream methylation \
  --method array \
  --betas /data/GSE59685_betas.csv \
  --metadata /data/meta_GSE59685.csv \
  --outdir results/methylation/
```

The beta matrix should have probe IDs (e.g. `cg00000029`) as rows and GSM
accession IDs as columns — the standard format for GEO supplementary files.

---

## Options

| Flag | Default | Description |
|------|---------|-------------|
| `--threads N` | 4 | CPU threads passed to each tool |
| `--explain` | on | Show brief educational explanations before each step |
| `--no-explain` | — | Suppress explanations (useful on repeated runs) |

---

## Output

Each pipeline command produces:

```
results/
├── <sample_name>/
│   ├── trimmed/          fastp output (trimmed FASTQs + QC report)
│   ├── aligned/          BAM files
│   ├── salmon/           (RNA-seq) Salmon quant.sf
│   ├── peaks/            (ATAC-seq) MACS2 narrowPeak files
│   └── methylation/      (WGBS) Bismark CX cytosine reports
└── obama_matrix.csv      Final output — ready for OBAMA
```

**`obama_matrix.csv`** has this structure:

| geo_accession | disease.state | GENE1 | GENE2 | … |
|---------------|---------------|-------|-------|---|
| GSM123456     | disease       | 14.2  | 0.3   | … |
| GSM123457     | control       | 2.1   | 8.7   | … |

> **Note on group labels:** the values in `disease.state` must be exactly `disease` and
> `control`. This is a hard requirement of the OBAMA pipeline's comparison code.
> Do not rename them even if your experiment does not involve a disease/control design.

---

## Reference Indices

Indices must be built once before running the pipeline. See `scripts/prepare_sample_data.sh`
for the exact commands. Typical locations on a shared server:

| Index | Path |
|-------|------|
| STAR (hg38) | `/data/upstream/indices/star_hg38` |
| Salmon (hg38) | `/data/upstream/indices/salmon_hg38` |
| Bowtie2 (hg38) | `/data/upstream/indices/bowtie2_hg38/hg38` |
| Bismark (hg38) | `/data/upstream/indices/bismark_hg38` |

---

## Running Tests

```bash
conda activate upstream
pip install pytest
pytest tests/
```

---

## Getting Help

```bash
upstream --help
upstream rnaseq --help
```
