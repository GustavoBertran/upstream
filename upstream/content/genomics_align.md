# Alignment — BWA-MEM

**BWA-MEM** maps each DNA read to its position on the reference genome. Unlike the
splice-aware aligners used for RNA (STAR), BWA-MEM expects reads to map contiguously
— correct for genomic DNA, which has no introns spliced out.

**Read groups (`@RG`)** are attached here and matter more than they look. The
`SM:` (sample) tag becomes the **sample column name in the VCF**, so every read from
sample `tumor_01` is tagged `SM:tumor_01` and the variant caller reports its
genotypes under that name. Without read groups, merging and per-sample bookkeeping
break.

Paired-end reads (R1 + R2) align more accurately than single-end — knowing the
distance between mates helps place reads in repetitive regions and improves indel
detection — but single-end works too.

Watch the mapping rate: for a human sample against the human reference, expect the
large majority of reads to map. A low rate usually means the wrong reference genome
or a contaminated library.
