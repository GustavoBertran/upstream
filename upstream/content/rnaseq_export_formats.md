## Counts Matrix + coldata — for DESeq2 / edgeR / limma-voom

Besides the OBAMA matrix (samples as rows), the pipeline can export the layout the
common R differential-expression tools expect: a **feature × sample** count matrix
plus a small **coldata** table describing each sample.

```
counts_matrix.csv          coldata.csv
gene,S1,S2,S3              sample,condition
GENE1,142,8,17            S1,disease
GENE2,0,1203,55           S2,control
...                        S3,control
```

Counts are **raw** (Salmon `NumReads`, STAR `ReadsPerGene`) — never TPM/FPKM —
because DESeq2, edgeR, and limma-voom model count data and do their own
normalization. Salmon counts are rounded to integers so `DESeqDataSetFromMatrix()`
accepts them; for a fully rigorous Salmon→DESeq2 path, import per-sample `quant.sf`
with `tximport` instead.

All three tools read these same two files — only the constructor differs:

```r
counts  <- as.matrix(read.csv("counts_matrix.csv", row.names = 1))
coldata <- read.csv("coldata.csv", row.names = 1)
coldata$condition <- factor(coldata$condition, levels = c("control", "disease"))
stopifnot(colnames(counts) == rownames(coldata))   # same order

# DESeq2
library(DESeq2)
dds <- DESeqDataSetFromMatrix(counts, coldata, design = ~ condition)
res <- results(DESeq(dds))

# edgeR
library(edgeR)
y   <- DGEList(counts, group = coldata$condition)
y   <- calcNormFactors(y[filterByExpr(y), , keep.lib.sizes = FALSE])
et  <- exactTest(estimateDisp(y, model.matrix(~ coldata$condition)))

# limma-voom
library(limma); library(edgeR)
v   <- voom(DGEList(counts), model.matrix(~ coldata$condition))
fit <- eBayes(lmFit(v, model.matrix(~ coldata$condition)))
```

`condition` uses the same `disease`/`control` labels as OBAMA. Set the factor
levels so `control` is the reference and log-fold-changes read disease-vs-control.
