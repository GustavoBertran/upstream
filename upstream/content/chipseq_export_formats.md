## Consensus-Peak Count Matrix — for DESeq2 / edgeR

Differential binding (is a region more bound in disease vs control?) is count-based,
so the default ChIP-seq output is a **peak × sample** raw-count matrix + a **coldata**
table — the same inputs DESeq2/edgeR use for RNA-seq.

```
counts_matrix.csv                  coldata.csv
peak,K27ac_tum,K27ac_norm          sample,condition
chr1:9800-10250,142,17             K27ac_tum,disease
chr3:1.2e6-1.21e6,0,1203           K27ac_norm,control
```

Two steps run for `--format matrix` (the default):

1. **Consensus peak set** — every ChIP sample's peaks are merged into one
   non-overlapping region set, so all samples are measured over the same features.
2. **Read counting** — `deeptools multiBamSummary` counts each sample's filtered
   reads in each consensus region.

```r
counts  <- as.matrix(read.csv("counts_matrix.csv", row.names = 1))
coldata <- read.csv("coldata.csv", row.names = 1)
coldata$condition <- factor(coldata$condition, levels = c("control", "disease"))
library(DESeq2)
dds <- DESeqDataSetFromMatrix(counts, coldata, design = ~ condition)
res <- results(DESeq(dds))   # differential binding per peak
```

**Why this is the default (not OBAMA).** OBAMA's interpretation modules (Gene
Ontology, Enrichr, STRING) are keyed on gene symbols. Peak coordinates map to none
of them, so `--format obama` for ChIP yields a peak-score matrix that OBAMA's
feature-selection (MCO) can read but its gene-based analyses cannot interpret. To
use OBAMA meaningfully on ChIP data, first annotate peaks to genes (e.g. ChIPseeker)
and build a gene-level matrix — a separate step not done here.
