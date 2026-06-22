## M-value Matrix + coldata — for limma

Besides the OBAMA matrix, the pipeline can export the layout limma expects for
differential methylation: a **CpG/probe × sample** matrix of **M-values** plus a
**coldata** table.

```
mvalues_matrix.csv          coldata.csv
cpg,S1,S2,S3               sample,condition
chr1:10468,2.31,-1.04,...  S1,disease
chr1:10470,-3.10,0.88,...  S2,control
...                         S3,control
```

**Why M-values, not beta/percent?** The M-value is `log2(beta / (1 - beta))`. Beta
values (0–1) are heteroscedastic — their variance is compressed near 0 and 1 — which
violates the constant-variance assumption of linear models. M-values are
approximately homoscedastic, so limma's linear modeling and empirical-Bayes
moderation behave correctly. Report effect sizes in beta for interpretability, but
**test on M-values** (Du et al., 2010).

- **Array (450K/EPIC):** M-value computed directly from beta, with beta clamped to
  `[1e-3, 1-1e-3]` so fully (un)methylated probes don't produce ±Inf.
- **WGBS (Bismark):** M-value from read counts as `log2((meth + 1) / (unmeth + 1))`;
  the pseudocount avoids ±Inf at low coverage. Only CpGs covered in every sample are
  kept (no imputation).

```r
M       <- as.matrix(read.csv("mvalues_matrix.csv", row.names = 1))
coldata <- read.csv("coldata.csv", row.names = 1)
coldata$condition <- factor(coldata$condition, levels = c("control", "disease"))
stopifnot(colnames(M) == rownames(coldata))

library(limma)
design <- model.matrix(~ condition, data = coldata)
fit    <- eBayes(lmFit(M, design))
top    <- topTable(fit, coef = "conditiondisease", number = Inf)
```

`condition` uses the same `disease`/`control` labels as OBAMA; `control` is the
reference level, so positive logFC means higher methylation in disease.
