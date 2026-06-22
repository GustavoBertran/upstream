## Consensus-Peak Count Matrix — for DESeq2 / edgeR

Differential accessibility is count-based, like RNA-seq differential expression — so
the inputs are the same shape: a **peak × sample** raw-count matrix + a **coldata**
table.

```
counts_matrix.csv                     coldata.csv
peak,S1,S2,S3                        sample,condition
chr1:9800-10250,142,8,17             S1,disease
chr1:28100-28640,0,1203,55           S2,control
...                                   S3,control
```

This is NOT the MACS2 score (the default OBAMA output). The MACS2 score is a
significance statistic, not a count, and can't go into a count model. Instead, two
steps run when you choose `--format matrix`:

1. **Consensus peak set** — every sample's peaks are merged into one non-overlapping
   set of regions, so all samples are measured over the *same* features.
2. **Read counting** — `deeptools multiBamSummary` counts the reads from each
   filtered BAM that fall in each consensus region (DiffBind-style).

```r
counts  <- as.matrix(read.csv("counts_matrix.csv", row.names = 1))
coldata <- read.csv("coldata.csv", row.names = 1)
coldata$condition <- factor(coldata$condition, levels = c("control", "disease"))
stopifnot(colnames(counts) == rownames(coldata))

library(DESeq2)
dds <- DESeqDataSetFromMatrix(counts, coldata, design = ~ condition)
res <- results(DESeq(dds))   # differential accessibility per peak
# edgeR works the same way (DGEList → estimateDisp → exactTest/glmQLFit)
```

Peak IDs are `chr:start-end` of the consensus region (also saved as
`consensus_peaks.bed`). `condition` uses the OBAMA `disease`/`control` labels, with
`control` as the reference level.
