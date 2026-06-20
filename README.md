# upstream

HTS preprocessing pipeline — trim, align, quantify, and produce
[OBAMA](https://github.com/AOG-Lab/OBAMA)-compatible output matrices.

See **[WELCOME.md](WELCOME.md)** for installation and usage instructions.

## Commands

```
upstream qc          FastQC + MultiQC on all samples
upstream rnaseq      QC → fastp → STAR → Salmon → OBAMA matrix
upstream atacseq     QC → fastp → Bowtie2 → filter → MACS2 → OBAMA matrix
upstream methylation QC → fastp → Bismark → methylation extract → OBAMA matrix
```

All commands accept `--no-explain` to suppress the step-by-step explanations shown by default.
