#!/usr/bin/env bash
# prepare_sample_data.sh
#
# Downloads and subsamples real GEO/SRA datasets for upstream tutorial use.
# Run as admin on the server.  Output goes to /data/upstream/shared/.
#
# Data source: GEO → SRA.  Find SRR accessions on a GEO sample page (GSM…)
# under "SRA Experiments", then fill in the variables below.
#
# Requirements: sra-tools (fasterq-dump, prefetch), seqtk
#   conda install -c bioconda sra-tools seqtk
#
# After running this script, students can generate their own samples.csv with:
#   upstream download --track <track> --outdir <workdir>

set -euo pipefail

SHARED=/data/upstream/shared
SEED=42
NREADS=1000000   # reads per sample; 1 M paired reads finishes in a ~1–3 h lab session

mkdir -p "${SHARED}/fastq/rnaseq"
mkdir -p "${SHARED}/fastq/atacseq"
mkdir -p "${SHARED}/fastq/methylation"
mkdir -p "${SHARED}/fastq/qc"

log() { echo "[$(date '+%H:%M:%S')] $*"; }

# ── SRR accessions ────────────────────────────────────────────────────────
# Fill these in from GEO before running.
# Each variable holds one SRR run accession (e.g. SRR1234567).
# For multi-run samples, concatenate FASTQ after fasterq-dump.

RNASEQ_DISEASE_SRR=""   # e.g. "SRR1234567"
RNASEQ_CONTROL_SRR=""   # e.g. "SRR1234568"

ATAC_DISEASE_SRR=""
ATAC_CONTROL_SRR=""

METHYL_DISEASE_SRR=""
METHYL_CONTROL_SRR=""

# ── Helpers ───────────────────────────────────────────────────────────────

check_srr() {
    local var_name="$1" srr="$2"
    if [[ -z "${srr}" ]]; then
        log "WARNING: ${var_name} is not set — skipping this sample."
        return 1
    fi
    return 0
}

fetch_and_subsample() {
    # Usage: fetch_and_subsample SRR_ACCESSION OUTDIR SAMPLE_NAME
    local srr="$1" outdir="$2" name="$3"
    local r1_final="${outdir}/${name}_R1.fastq.gz"
    local r2_final="${outdir}/${name}_R2.fastq.gz"

    if [[ -f "${r1_final}" && -f "${r2_final}" ]]; then
        log "Already have ${name} FASTQ — skipping."
        return
    fi

    log "Fetching ${srr} (${name}) via fasterq-dump …"
    fasterq-dump \
        --split-files \
        --outdir "${outdir}" \
        --progress \
        "${srr}"

    local raw_r1="${outdir}/${srr}_1.fastq"
    local raw_r2="${outdir}/${srr}_2.fastq"

    if [[ ! -f "${raw_r1}" || ! -f "${raw_r2}" ]]; then
        log "ERROR: fasterq-dump did not produce expected files for ${srr}"
        exit 1
    fi

    log "Subsampling ${name} to ${NREADS} reads …"
    seqtk sample -s "${SEED}" "${raw_r1}" "${NREADS}" | gzip > "${r1_final}"
    seqtk sample -s "${SEED}" "${raw_r2}" "${NREADS}" | gzip > "${r2_final}"

    rm -f "${raw_r1}" "${raw_r2}"
    log "Done → ${r1_final##*/}, ${r2_final##*/}"
}

# ── RNA-seq ───────────────────────────────────────────────────────────────

log "=== RNA-seq ==="
if check_srr "RNASEQ_DISEASE_SRR" "${RNASEQ_DISEASE_SRR}"; then
    fetch_and_subsample "${RNASEQ_DISEASE_SRR}" "${SHARED}/fastq/rnaseq" "disease"
fi
if check_srr "RNASEQ_CONTROL_SRR" "${RNASEQ_CONTROL_SRR}"; then
    fetch_and_subsample "${RNASEQ_CONTROL_SRR}" "${SHARED}/fastq/rnaseq" "control"
fi

# QC module reuses rnaseq FASTQ
if [[ -f "${SHARED}/fastq/rnaseq/disease_R1.fastq.gz" ]]; then
    ln -sf "${SHARED}/fastq/rnaseq/disease_R1.fastq.gz" "${SHARED}/fastq/qc/disease_R1.fastq.gz"
    ln -sf "${SHARED}/fastq/rnaseq/disease_R2.fastq.gz" "${SHARED}/fastq/qc/disease_R2.fastq.gz"
    ln -sf "${SHARED}/fastq/rnaseq/control_R1.fastq.gz" "${SHARED}/fastq/qc/control_R1.fastq.gz"
    ln -sf "${SHARED}/fastq/rnaseq/control_R2.fastq.gz" "${SHARED}/fastq/qc/control_R2.fastq.gz"
fi

# ── ATAC-seq ──────────────────────────────────────────────────────────────

log "=== ATAC-seq ==="
if check_srr "ATAC_DISEASE_SRR" "${ATAC_DISEASE_SRR}"; then
    fetch_and_subsample "${ATAC_DISEASE_SRR}" "${SHARED}/fastq/atacseq" "disease"
fi
if check_srr "ATAC_CONTROL_SRR" "${ATAC_CONTROL_SRR}"; then
    fetch_and_subsample "${ATAC_CONTROL_SRR}" "${SHARED}/fastq/atacseq" "control"
fi

# ── Methylation ───────────────────────────────────────────────────────────

log "=== Methylation ==="
if check_srr "METHYL_DISEASE_SRR" "${METHYL_DISEASE_SRR}"; then
    fetch_and_subsample "${METHYL_DISEASE_SRR}" "${SHARED}/fastq/methylation" "disease"
fi
if check_srr "METHYL_CONTROL_SRR" "${METHYL_CONTROL_SRR}"; then
    fetch_and_subsample "${METHYL_CONTROL_SRR}" "${SHARED}/fastq/methylation" "control"
fi

# ── Reference indices (reminder — build separately) ───────────────────────

cat << 'EOF'

=======================================================================
NEXT STEP: Build reference indices
=======================================================================

Store indices at /data/upstream/shared/indices/ and set permissions:

1. STAR index (hg38, ~30 GB RAM):
   STAR --runMode genomeGenerate \
        --genomeDir /data/upstream/shared/indices/star_hg38 \
        --genomeFastaFiles /path/to/hg38.fa \
        --sjdbGTFfile /path/to/gencode.v44.annotation.gtf \
        --runThreadN 8

2. Salmon index:
   salmon index \
        --transcripts /path/to/gencode.v44.transcripts.fa.gz \
        --index /data/upstream/shared/indices/salmon_hg38 \
        --gencode --threads 8

3. Bowtie2 index (hg38):
   bowtie2-build /path/to/hg38.fa \
        /data/upstream/shared/indices/bowtie2_hg38/hg38

4. Bismark genome (hg38, ~100 GB disk):
   bismark_genome_preparation /path/to/hg38_dir/
   # Point --bismark-genome to: /data/upstream/shared/indices/bismark_hg38

5. Set permissions:
   chmod -R 755 /data/upstream/shared/
   chmod -R 775 /data/upstream/students/

EOF

log "Sample data preparation complete."
