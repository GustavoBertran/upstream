## Alignment — Bowtie2

Trimmed reads are aligned to the reference genome with Bowtie2 (`--very-sensitive`)
and the output is coordinate-sorted and indexed with samtools. ChIP-seq needs only
a standard genome alignment — there is no splice awareness (unlike RNA-seq) and no
Tn5 adjustment (unlike ATAC-seq).

Both the ChIP (IP) samples and any input/control samples are aligned the same way;
the input simply provides the background read distribution that MACS2 subtracts.
