# Variant Calling — bcftools or GATK

This step finds where each sample differs from the reference and writes a **VCF** —
one line per variant: chromosome, position, reference allele, alternate allele,
quality, and per-sample genotype. Two callers are available (`--caller`):

**bcftools (default, lightweight).** **bcftools mpileup** walks the genome and, at
each position, summarizes the aligned bases (the "pileup"); **bcftools call -mv**
then calls **variants only** (`-v`) with the multiallelic model (`-m`). Fast, few
dependencies — ideal for learning.

**GATK HaplotypeCaller (`--caller gatk`, the field standard).** Rather than reading
the pileup column by column, HaplotypeCaller **locally re-assembles** each active
region from the reads and realigns them, which improves calls around indels and in
messy regions. It is the **"lingua franca"** of germline calling — what most
published pipelines and clinical workflows use. It needs the reference `.fai` **and
a sequence dictionary** (`.dict`; `samtools dict ref.fa -o ref.dict`), which
`index-help --tool bwa` sets up.

Both callers here run **per sample**, then (optionally) merge the per-sample VCFs
into one `cohort.vcf.gz`. That is a teaching simplification of **joint genotyping**.

> **Gold standard for production:** the full **GATK Best Practices** pipeline — base
> quality score recalibration (BQSR), then HaplotypeCaller in **GVCF** mode per
> sample, then **GenotypeGVCFs** to genotype the cohort jointly, then filtering
> (VQSR or hard filters). That distinguishes a confident `0/0` from a `./.` "not
> looked at" — which merging independently-called VCFs cannot. `--caller gatk` uses
> the same caller; the GVCF/joint steps are the remaining production layer.
