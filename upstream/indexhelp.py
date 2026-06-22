"""Recommended commands for obtaining/building a reference index.

This module only *prints guidance* — it never downloads or builds anything (index
builds are large and machine-specific: STAR needs ~30 GB RAM, Bismark ~100 GB disk
for human). Reference URLs use GENCODE; bump the release if you want a newer one.
"""
from __future__ import annotations

TOOLS = ("salmon", "star", "bowtie2", "bismark")

GENOMES = {
    "human": {
        "label": "human GRCh38 (GENCODE release 44)",
        "base": "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44",
        "genome": "GRCh38.primary_assembly.genome.fa.gz",
        "transcripts": "gencode.v44.transcripts.fa.gz",
        "gtf": "gencode.v44.primary_assembly.annotation.gtf.gz",
    },
    "mouse": {
        "label": "mouse GRCm39 (GENCODE release M34)",
        "base": "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M34",
        "genome": "GRCm39.primary_assembly.genome.fa.gz",
        "transcripts": "gencode.vM34.pc_transcripts.fa.gz",
        "gtf": "gencode.vM34.primary_assembly.annotation.gtf.gz",
    },
}


def recommend(tool: str, genome: str = "human") -> str:
    g = GENOMES[genome]
    base, fa_gz, tx_gz, gtf_gz = g["base"], g["genome"], g["transcripts"], g["gtf"]
    fa, gtf = fa_gz[:-3], gtf_gz[:-3]  # uncompressed names

    if tool == "salmon":
        return (
            f"# Salmon decoy-aware index — {g['label']} (~4 GB RAM, ~45 min)\n"
            f"curl -O {base}/{tx_gz}\n"
            f"curl -O {base}/{fa_gz}\n"
            f"gzip -dc {fa_gz} | grep '^>' | cut -d ' ' -f1 | sed 's/>//' > decoys.txt\n"
            f"cat {tx_gz} {fa_gz} > gentrome.fa.gz\n"
            f"salmon index -t gentrome.fa.gz -d decoys.txt -i salmon_index --gencode -k 31 -p 8\n"
            f"# then:  upstream rnaseq --aligner salmon --salmon-index ./salmon_index ..."
        )
    if tool == "star":
        return (
            f"# STAR index — {g['label']}  (NEEDS ~30 GB RAM to build; use a server)\n"
            f"curl -O {base}/{fa_gz}\n"
            f"curl -O {base}/{gtf_gz}\n"
            f"gzip -d {fa_gz} {gtf_gz}\n"
            f"mkdir -p star_index\n"
            f"STAR --runMode genomeGenerate --genomeDir star_index \\\n"
            f"     --genomeFastaFiles {fa} --sjdbGTFfile {gtf} --runThreadN 8\n"
            f"# then:  upstream rnaseq --aligner star --star-index ./star_index ..."
        )
    if tool == "bowtie2":
        return (
            f"# Bowtie2 index — {g['label']}  (a few GB RAM)\n"
            f"curl -O {base}/{fa_gz}\n"
            f"gzip -d {fa_gz}\n"
            f"mkdir -p bowtie2_index\n"
            f"bowtie2-build --threads 8 {fa} bowtie2_index/genome\n"
            f"# then:  upstream atacseq --bowtie2-index ./bowtie2_index/genome ..."
        )
    if tool == "bismark":
        return (
            f"# Bismark genome — {g['label']}  (~100 GB disk for human; use a server)\n"
            f"mkdir -p bismark_genome && cd bismark_genome\n"
            f"curl -O {base}/{fa_gz}\n"
            f"gzip -d {fa_gz}            # the FASTA must sit inside this folder\n"
            f"bismark_genome_preparation .\n"
            f"# then:  upstream methylation --method wgbs --bismark-genome ./bismark_genome ..."
        )
    raise ValueError(f"unknown tool: {tool}")
