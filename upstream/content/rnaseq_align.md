# Genome Alignment — STAR (RNA-seq)

mRNA is spliced: introns are removed and exons are joined. When you sequence a cDNA library,
some reads span exon–exon junctions that are far apart in the genome. Standard aligners fail
on these reads. STAR handles them by detecting splice junctions from the data itself.

**What to watch for in `_Log.final.out`:**
- Uniquely mapped reads > 70% — healthy
- < 60% — usually means wrong reference genome, degraded RNA, or heavy contamination
- High % multi-mapped — can indicate repetitive sequence or gene family members; normal < 20%

The output BAM is coordinate-sorted, which is required for Salmon, IGV viewing, and most
downstream tools.
