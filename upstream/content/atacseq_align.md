# Genome Alignment — Bowtie2 (ATAC-seq)

ATAC-seq reads come from genomic DNA, not spliced mRNA — no exon–intron junctions, so
STAR's splice-aware alignment is unnecessary. Bowtie2 is faster and designed for this case.

**Key flags:**
- `--no-mixed / --no-discordant` — discard pairs where only one mate aligns or mates
  align in unexpected orientations (these are likely artefacts)
- `-X 2000` — allow fragment sizes up to 2 kb to capture nucleosome-sized fragments

The Bowtie2 output is piped directly to `samtools sort`, avoiding a large intermediate SAM
file. The BAM is then indexed for all downstream filtering and peak-calling steps.

Expect > 70% overall alignment rate for a good ATAC-seq library.
