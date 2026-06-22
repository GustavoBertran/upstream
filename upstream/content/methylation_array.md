## Illumina 450K / EPIC Array — OBAMA Matrix

This step reads a beta-value matrix downloaded from GEO and reformats it into an OBAMA-compatible CSV.

**Beta values** (range 0–1) represent the proportion of methylated DNA at each CpG locus: 0 = fully unmethylated, 1 = fully methylated. OBAMA treats these as continuous features and compares disease vs. control groups across all probes.

**Probe IDs** (e.g., `cg06432309`) are Illumina-specific identifiers for ~485,000 CpGs on the 450K array (~850,000 on EPIC). After MCO analysis, annotate probe IDs back to genomic coordinates and gene context using `minfi` in R:

```r
library(minfi)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
anno <- as.data.frame(getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19))
```

**NA filtering**: Probes with missing values in any selected sample are dropped. Imputing 0.0 for failed probes would create a false "fully unmethylated" signal in OBAMA — the same reason WGBS only uses CpGs covered in every sample.
