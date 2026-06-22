## Alignment Filtering

Before peak calling, alignments are cleaned up:

- **Mitochondrial / non-standard contigs removed** — reads are restricted to the
  autosomes plus X/Y. Mitochondrial and unplaced-contig reads are background noise.
- **Duplicates removed** (`samtools markdup -r`) — PCR/optical duplicates would
  create artificial pileups that look like enrichment.
- **MAPQ ≥ 30** — keeps confidently, uniquely mapped reads; multi-mappers in
  repetitive regions produce false peaks.

The same filtering is applied to input/control samples so the IP and its background
are compared on equal footing.
