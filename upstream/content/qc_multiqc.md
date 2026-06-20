# MultiQC — Aggregate QC Reports

MultiQC scans the FastQC output directory and combines everything into a single interactive
HTML report. One heatmap lets you spot outlier samples instantly — a row that is grey while
all others are green tells you exactly where to investigate.

**General statistics table** — read counts, % duplicates, average quality, % GC per sample.
Any sample that is an outlier in any column warrants a closer look.

**Per-base quality heatmap** — compact view of quality across all samples simultaneously.

**Adapter content bars** — height shows what fraction of reads carry adapter contamination.
