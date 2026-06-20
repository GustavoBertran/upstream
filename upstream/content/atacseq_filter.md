# Filtering — Mitochondrial Reads, Duplicates, MAPQ (ATAC-seq)

Three filtering steps are run in sequence:

**1. Remove mitochondrial reads (chrM)**
Mitochondria have no nucleosomes. Tn5 cuts them freely, so 30–80% of ATAC-seq reads
can come from chrM. These reads tell you nothing about chromatin accessibility and would
dominate peak calling. We keep only nuclear chromosomes (chr1–chr22, chrX, chrY).

**2. Remove PCR duplicates**
PCR amplification creates identical copies of the same original fragment. These are not
independent observations — counting them inflates accessibility at those positions.
`samtools markdup -r` identifies read pairs with identical start/end coordinates and
removes all but one.

**3. Filter MAPQ ≥ 30**
MAPQ (mapping quality) reflects aligner confidence in a read's placement. Multi-mapping
reads (equally plausible at multiple positions) get low MAPQ. Keeping MAPQ ≥ 30 means
retaining reads the aligner is ≥ 99.9% confident are correctly placed.
