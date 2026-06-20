#!/usr/bin/env bash
# prepare_sample_data.sh
#
# Downloads and subsamples real public HTS datasets for upstream tutorial use.
# Run as admin on the server. Output goes to /data/upstream/shared/.
#
# Datasets chosen:
#   RNA-seq:     MCF-7 breast cancer cell line, poly-A RNA (ENCODE)
#                  disease: ENCFF000EFG (MCF-7)   → simulates "disease"
#                  control: ENCFF000EFH (MCF-7 parental) → simulates "control"
#   ATAC-seq:    MCF-7 ATAC-seq (ENCODE ENCSR000DZL)
#                  disease: ENCFF828ZPN
#                  control: ENCFF585IFG
#   Methylation: MCF-7 WGBS (ENCODE ENCSR765JPC)
#                  disease: ENCFF721JMB
#                  control: ENCFF271RWW
#
# All files are subsampled to ~1,000,000 paired reads using seqtk so that
# each track completes in 20–60 minutes on modest shared compute.
#
# Requirements: seqtk, wget or curl, samtools (for ATAC/methylation ENCODE files)
# The ENCODE files are BAMs; we convert to FASTQ with samtools fastq.

set -euo pipefail

SHARED=/data/upstream/shared
SEED=42
NREADS=1000000   # per sample; 1M paired reads ~= 300 MB gzipped FASTQ

mkdir -p "${SHARED}/fastq/shared_qc"
mkdir -p "${SHARED}/fastq/rnaseq"
mkdir -p "${SHARED}/fastq/atacseq"
mkdir -p "${SHARED}/fastq/methylation"
mkdir -p "${SHARED}/indices"

log() { echo "[$(date '+%H:%M:%S')] $*"; }

# ── Helpers ──────────────────────────────────────────────────────────────

subsample_fastq_pair() {
    local r1_in="$1" r2_in="$2" r1_out="$3" r2_out="$4"
    log "Subsampling $(basename "${r1_in}") to ${NREADS} reads …"
    seqtk sample -s "${SEED}" "${r1_in}" "${NREADS}" | gzip > "${r1_out}"
    seqtk sample -s "${SEED}" "${r2_in}" "${NREADS}" | gzip > "${r2_out}"
    log "Done → $(basename "${r1_out}"), $(basename "${r2_out}")"
}

bam_to_fastq_pair() {
    local bam="$1" r1_out="$2" r2_out="$3"
    local tmp_r1="${bam%.bam}_R1.fastq"
    local tmp_r2="${bam%.bam}_R2.fastq"
    log "Converting BAM to FASTQ: $(basename "${bam}") …"
    samtools sort -n "${bam}" -o "${bam%.bam}.namesorted.bam"
    samtools fastq -1 "${tmp_r1}" -2 "${tmp_r2}" \
        -0 /dev/null -s /dev/null -n \
        "${bam%.bam}.namesorted.bam"
    gzip "${tmp_r1}" "${tmp_r2}"
    echo "${tmp_r1}.gz" "${tmp_r2}.gz"
}

fetch() {
    local url="$1" dest="$2"
    if [[ -f "${dest}" ]]; then
        log "Already have $(basename "${dest}"), skipping download."
        return
    fi
    log "Downloading $(basename "${dest}") …"
    wget --quiet --show-progress -O "${dest}" "${url}" || \
    curl -L --progress-bar -o "${dest}" "${url}"
}

# ── RNA-seq ──────────────────────────────────────────────────────────────
# ENCODE accession: ENCSR000CRB (MCF-7 polyA RNA-seq, paired-end)
# Replace these URLs with accessions confirmed appropriate for the course organism.

log "=== RNA-seq sample data ==="

RNA_DISEASE_URL="https://www.encodeproject.org/files/ENCFF000EFG/@@download/ENCFF000EFG.fastq.gz"
RNA_CONTROL_URL="https://www.encodeproject.org/files/ENCFF000EFH/@@download/ENCFF000EFH.fastq.gz"

# NOTE: ENCODE RNA-seq files are typically single-end at older accessions.
# For paired-end, you may need two separate accessions per sample.
# Adjust the URLs and subsample commands below once you confirm the accessions.

RNA_DISEASE_R1_RAW="${SHARED}/fastq/rnaseq/.raw_disease_R1.fastq.gz"
RNA_DISEASE_R2_RAW="${SHARED}/fastq/rnaseq/.raw_disease_R2.fastq.gz"
RNA_CONTROL_R1_RAW="${SHARED}/fastq/rnaseq/.raw_control_R1.fastq.gz"
RNA_CONTROL_R2_RAW="${SHARED}/fastq/rnaseq/.raw_control_R2.fastq.gz"

# fetch "${RNA_DISEASE_R1_URL}" "${RNA_DISEASE_R1_RAW}"   # uncomment with real URLs
# fetch "${RNA_DISEASE_R2_URL}" "${RNA_DISEASE_R2_RAW}"

if [[ -f "${RNA_DISEASE_R1_RAW}" && -f "${RNA_DISEASE_R2_RAW}" ]]; then
    subsample_fastq_pair \
        "${RNA_DISEASE_R1_RAW}" "${RNA_DISEASE_R2_RAW}" \
        "${SHARED}/fastq/rnaseq/disease_R1.fastq.gz" \
        "${SHARED}/fastq/rnaseq/disease_R2.fastq.gz"
    subsample_fastq_pair \
        "${RNA_CONTROL_R1_RAW}" "${RNA_CONTROL_R2_RAW}" \
        "${SHARED}/fastq/rnaseq/control_R1.fastq.gz" \
        "${SHARED}/fastq/rnaseq/control_R2.fastq.gz"
    # Copy to shared_qc/ (Module 0 uses whatever track's data is available)
    ln -sf "${SHARED}/fastq/rnaseq/disease_R1.fastq.gz" "${SHARED}/fastq/shared_qc/disease_R1.fastq.gz"
    ln -sf "${SHARED}/fastq/rnaseq/disease_R2.fastq.gz" "${SHARED}/fastq/shared_qc/disease_R2.fastq.gz"
    ln -sf "${SHARED}/fastq/rnaseq/control_R1.fastq.gz" "${SHARED}/fastq/shared_qc/control_R1.fastq.gz"
    ln -sf "${SHARED}/fastq/rnaseq/control_R2.fastq.gz" "${SHARED}/fastq/shared_qc/control_R2.fastq.gz"
else
    log "WARNING: RNA-seq raw files not found. Download them and re-run."
fi

# ── ATAC-seq ─────────────────────────────────────────────────────────────
# ENCODE accession: ENCSR000DZL (MCF-7 ATAC-seq)

log "=== ATAC-seq sample data ==="

ATAC_DISEASE_BAM_URL="https://www.encodeproject.org/files/ENCFF828ZPN/@@download/ENCFF828ZPN.bam"
ATAC_CONTROL_BAM_URL="https://www.encodeproject.org/files/ENCFF585IFG/@@download/ENCFF585IFG.bam"

ATAC_DISEASE_BAM="${SHARED}/fastq/atacseq/.raw_disease.bam"
ATAC_CONTROL_BAM="${SHARED}/fastq/atacseq/.raw_control.bam"

fetch "${ATAC_DISEASE_BAM_URL}" "${ATAC_DISEASE_BAM}"
fetch "${ATAC_CONTROL_BAM_URL}" "${ATAC_CONTROL_BAM}"

for SAMPLE in disease control; do
    BAM="${SHARED}/fastq/atacseq/.raw_${SAMPLE}.bam"
    R1="${SHARED}/fastq/atacseq/${SAMPLE}_R1.fastq.gz"
    R2="${SHARED}/fastq/atacseq/${SAMPLE}_R2.fastq.gz"
    if [[ -f "${BAM}" && ! -f "${R1}" ]]; then
        read -r TMP_R1 TMP_R2 <<< "$(bam_to_fastq_pair "${BAM}" "${R1}" "${R2}")"
        subsample_fastq_pair "${TMP_R1}" "${TMP_R2}" "${R1}" "${R2}"
        rm -f "${TMP_R1}" "${TMP_R2}"
    fi
done

# ── Methylation ──────────────────────────────────────────────────────────
# ENCODE accession: ENCSR765JPC (MCF-7 WGBS)

log "=== Methylation (WGBS) sample data ==="

WGBS_DISEASE_BAM_URL="https://www.encodeproject.org/files/ENCFF721JMB/@@download/ENCFF721JMB.bam"
WGBS_CONTROL_BAM_URL="https://www.encodeproject.org/files/ENCFF271RWW/@@download/ENCFF271RWW.bam"

WGBS_DISEASE_BAM="${SHARED}/fastq/methylation/.raw_disease.bam"
WGBS_CONTROL_BAM="${SHARED}/fastq/methylation/.raw_control.bam"

fetch "${WGBS_DISEASE_BAM_URL}" "${WGBS_DISEASE_BAM}"
fetch "${WGBS_CONTROL_BAM_URL}" "${WGBS_CONTROL_BAM}"

for SAMPLE in disease control; do
    BAM="${SHARED}/fastq/methylation/.raw_${SAMPLE}.bam"
    R1="${SHARED}/fastq/methylation/${SAMPLE}_R1.fastq.gz"
    R2="${SHARED}/fastq/methylation/${SAMPLE}_R2.fastq.gz"
    if [[ -f "${BAM}" && ! -f "${R1}" ]]; then
        read -r TMP_R1 TMP_R2 <<< "$(bam_to_fastq_pair "${BAM}" "${R1}" "${R2}")"
        subsample_fastq_pair "${TMP_R1}" "${TMP_R2}" "${R1}" "${R2}"
        rm -f "${TMP_R1}" "${TMP_R2}"
    fi
done

# ── Reference indices (reminder — these must be built separately) ─────────

cat << 'EOF'

=======================================================================
NEXT STEP: Build reference indices
=======================================================================

The following indices must be built once by the admin and stored at:
  /data/upstream/shared/indices/

1. STAR index (hg38, requires ~30 GB RAM):
   STAR --runMode genomeGenerate \
        --genomeDir /data/upstream/shared/indices/star_hg38 \
        --genomeFastaFiles /path/to/hg38.fa \
        --sjdbGTFfile /path/to/gencode.v44.annotation.gtf \
        --runThreadN 8

2. Salmon index (from GENCODE transcriptome FASTA):
   salmon index \
        --transcripts /path/to/gencode.v44.transcripts.fa.gz \
        --index /data/upstream/shared/indices/salmon_hg38 \
        --gencode --threads 8

3. Bowtie2 index (hg38):
   bowtie2-build /path/to/hg38.fa \
        /data/upstream/shared/indices/bowtie2_hg38/hg38

4. Bismark index (hg38, requires ~100 GB disk for C-T conversion):
   bismark_genome_preparation /path/to/hg38_dir/ \
        --genomic_composition
   # Point --genome in bismark commands to: /data/upstream/shared/indices/bismark_hg38

5. Set permissions so all students can read but not write:
   chmod -R 755 /data/upstream/shared/
   chmod -R 775 /data/upstream/students/

EOF

log "Sample data preparation complete."
