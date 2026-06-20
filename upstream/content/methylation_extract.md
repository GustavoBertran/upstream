# Methylation Extraction — bismark_methylation_extractor

The Bismark BAM contains aligned reads. This step reads each read and, for every cytosine
position in the genome, records whether that read shows C (methylated) or T (unconverted,
unmethylated). The `--cytosine_report` flag produces a genome-wide table of every covered
cytosine with its coverage depth and methylation count.

**CpG vs CHH vs CHG:**
- **CpG** — the main target of DNA methylation in mammals. Methylated CpG at promoters
  silences gene expression; enhancers gain methylation when inactive.
- **CHH / CHG** — nearly always unmethylated in mammals (H = A, T, or C). Very high
  non-CpG methylation (> 2%) suggests incomplete bisulfite conversion or plant DNA
  contamination.

`--CX_context` reports all three contexts so you can use non-CpG methylation as an
internal quality check.

The CX report is used in the final step to build the OBAMA matrix: one column per CpG
position, one row per sample, values are percent methylation (0–100).
