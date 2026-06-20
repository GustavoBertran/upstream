# FastQC — Raw Read Quality Check

FastQC scans each FASTQ file and reports per-base quality scores, adapter contamination,
GC content, and duplication levels. Catching problems here saves hours of wasted analysis.

**Per-base Phred quality** — Q30 = 1 in 1,000 error chance; Q20 = 1 in 100. A gradual drop
toward read ends is normal. A crash to Q10 at position 20 is not.

**Adapter content** — if your insert is shorter than the read length, the sequencer reads into
the adapter. Aligners will not map adapter sequence to the genome; trimming (next step) fixes it.

**GC content** — human reads should cluster around 41%. A bimodal distribution often means
contamination with another organism.

**Duplication** — context matters: high duplication in RNA-seq can indicate low-complexity
libraries; in ATAC-seq it is expected and removed downstream.

Open the HTML reports when done. Green ✓ = fine, yellow ⚠ = check it, red ✗ = problem.
