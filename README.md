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

**RNA-seq:**
```bash
upstream rnaseq \
  --samples samples.csv \
  --star-index   /path/to/star_hg38 \
  --salmon-index /path/to/salmon_hg38 \
  --outdir results/rnaseq/
```

**ATAC-seq:**
```bash
upstream atacseq \
  --samples samples.csv \
  --bowtie2-index /path/to/bowtie2/hg38 \
  --outdir results/atacseq/
```

**DNA Methylation (WGBS/RRBS):**
```bash
upstream methylation \
  --samples samples.csv \
  --bismark-genome /path/to/bismark_hg38 \
  --outdir results/methylation/
```

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
