# Transcript Quantification — Salmon

Salmon uses quasi-mapping (not full alignment) to estimate expression levels. It handles
multi-mapping reads with an EM algorithm, giving more accurate isoform estimates than
simple alignment-and-count approaches.

**TPM (Transcripts Per Million)** — the output unit:
- Corrects for transcript length (longer genes capture more reads at equal expression)
- Normalises for sequencing depth
- Sum = exactly 1,000,000 per sample, making samples directly comparable

Watch the mapping rate printed to the terminal: > 75% is excellent, < 50% suggests the
wrong transcriptome index or heavy contamination.

The `quant.sf` output feeds directly into DESeq2/edgeR via the tximport/tximeta R packages.
