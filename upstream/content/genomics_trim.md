# Read Trimming — fastp

Before aligning DNA reads to the genome, **fastp** removes sequencing adapters and
low-quality bases. Leftover adapter sequence and noisy read ends create false
mismatches at alignment time, which downstream looks like spurious variants.

This is the same trimming step every track starts with — for variant calling its
job is simply to give BWA-MEM clean reads so that the mismatches that *do* survive
are real differences from the reference, not artifacts.

Watch the percentage of reads passing the filter: a healthy DNA library keeps the
large majority of its reads.
