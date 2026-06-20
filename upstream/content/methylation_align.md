# Bisulfite Alignment — Bismark

Standard aligners cannot handle bisulfite reads: after conversion, every unmethylated C
looks like a T, so the read would appear to have thousands of mismatches against the genome.

Bismark solves this by:
1. **Converting the reference** in silico — creating a C→T version and a G→A version of
   every chromosome
2. **Converting the reads** the same way
3. **Aligning converted reads** to converted genomes using Bowtie2 under the hood
4. **Reporting positions** in original genome coordinates

This means every genomic C is treated as a wildcard that matches either C (methylated) or
T (unconverted = unmethylated) in the read.

**Bisulfite conversion efficiency** — reported in `_bismark_bt2_PE_report.txt`:
non-CpG cytosines in mammals are nearly always unmethylated, so they should all appear
as T. > 98% = excellent; < 95% = possible incomplete conversion, results may be unreliable.

Bismark is slower than most aligners (it runs Bowtie2 twice). Expect 10–30 min on
subsampled data with `-p 2`.
