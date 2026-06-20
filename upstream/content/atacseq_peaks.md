# Peak Calling — MACS2 (ATAC-seq)

A "peak" is a genomic region with significantly more reads than expected by chance —
where chromatin is open and Tn5 could insert. Peak callers identify these pileups
and return them as genomic intervals (BED-like format).

**The ATAC-seq shift:**
Raw reads align to where Tn5 cut. But Tn5 inserts at an offset from the cut site
(+4 bp on the + strand, −5 bp on the − strand). `--shift -100 --extsize 200` re-centres
reads around the actual Tn5 insertion site so peak summits mark the true centre of
open chromatin — important for transcription factor footprinting downstream.

`--nomodel` disables MACS2's ChIP-seq fragment-size model, which would give wrong results
here. `--nomodel --shift -100 --extsize 200` is the standard ATAC-seq setting.

A good ATAC-seq experiment yields > 30,000 peaks (on full-depth data; subsampled data
will yield fewer). Fraction of reads in peaks (FRiP) > 20% indicates a clean library.
