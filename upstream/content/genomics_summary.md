# Variant Summary

This step indexes and (optionally) merges the per-sample VCFs and writes
`variant_summary.csv` — one row per sample with total variant records, SNP and indel
counts, and the **ts/tv ratio** (transitions ÷ transversions).

**ts/tv** is a quick quality signal: real human germline variants sit around
**2.0–2.1** genome-wide (higher, ~3, in exomes). A ratio near 0.5 — the value you'd
get from random base changes — suggests many false positives.

**Reading a VCF line:** `CHROM POS ID REF ALT QUAL FILTER INFO FORMAT <sample…>`.
A genotype of `0/1` is heterozygous (one reference, one alternate allele), `1/1` is
homozygous alternate.

> **A caveat about `cohort.vcf.gz`.** Because samples were called *independently* and
> then merged, a site where one sample had no call shows `./.` (missing) for that
> sample — **not** a confident `0/0` (homozygous reference). Only true joint
> genotyping (one pileup over all samples, or GATK GenotypeGVCFs) can distinguish
> "reference here" from "we didn't look here." Don't over-interpret missing
> genotypes in the merged cohort VCF.

This track stops at raw SNVs and indels — it does not recalibrate, filter, annotate
(VEP/SnpEff), or call structural/copy-number variants.
