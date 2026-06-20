# Adapter Trimming — fastp (ATAC-seq)

ATAC-seq uses the Tn5 transposase to insert **Nextera adapters** directly into open chromatin.
Short insert sizes are common: nucleosome-free regions produce sub-150 bp fragments that are
almost entirely adapter sequence in a 150 bp read. These short reads are valuable — they
come from the most accessible chromatin — so we trim adapters rather than discard the reads.

High adapter contamination (> 30%) is **normal and expected** for ATAC-seq. After trimming,
you should see a multimodal read-length distribution reflecting nucleosome-free (~180 bp)
and mono-nucleosome (~360 bp) fragments.
