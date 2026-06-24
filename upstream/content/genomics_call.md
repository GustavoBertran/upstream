# Variant Calling — bcftools

**bcftools mpileup** walks the genome and, at each position, summarizes the aligned
bases (the "pileup"); **bcftools call -mv** then decides where the sample differs
from the reference, emitting **variants only** (`-v`) using the multiallelic model
(`-m`). The result is a **VCF** — one line per variant: chromosome, position,
reference allele, alternate allele, quality, and per-sample genotype.

This track calls each sample **independently**, then (optionally) merges the
per-sample VCFs into one `cohort.vcf.gz`. That is the most common scriptable
approach and is ideal for learning, but it is a simplification of **joint
genotyping**.

> **Gold standard for production:** the **GATK Best Practices** pipeline — base
> quality score recalibration (BQSR), then **HaplotypeCaller** in GVCF mode per
> sample, then **GenotypeGVCFs** to genotype the cohort jointly, then variant
> filtering (VQSR or hard filters). It is heavier (needs known-sites resources and
> more RAM/time), which is why this teaching track uses bcftools instead. The
> concepts you learn here transfer directly.
