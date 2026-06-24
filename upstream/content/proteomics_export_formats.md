## Proteomics Export — limma matrix

This track writes two files for differential-abundance analysis in R with **limma**
(the standard tool for proteomics, just as DESeq2/edgeR are for RNA-seq counts):

- **`proteins_matrix.csv`** — proteins as rows, samples as columns. Values are
  **normalized log2 intensities**. Blank cells are missing values (read as `NA` in R)
  — limma tolerates them.
- **`coldata.csv`** — `sample,condition`, where condition is `disease` or `control`.
  The sample names match the matrix columns exactly.

Load and test in R:

```r
library(limma)
mat <- as.matrix(read.csv("proteins_matrix.csv", row.names = 1))
col <- read.csv("coldata.csv")
design <- model.matrix(~ condition, data = col)   # control = reference level
fit <- eBayes(lmFit(mat, design))
topTable(fit, coef = 2)                            # disease vs control
```

limma models the log2 intensities directly with an empirical-Bayes moderated
t-test, which is well-suited to the small sample sizes and many features typical of
proteomics. If you skipped normalization here (`--normalize none`), do it in R with
`normalizeBetweenArrays(mat)` before fitting.
