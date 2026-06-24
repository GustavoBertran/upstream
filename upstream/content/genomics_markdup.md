# Mark Duplicates — samtools markdup

**PCR duplicates** are multiple read pairs that came from the *same* original DNA
fragment, copied during library amplification. They are not independent evidence —
if you count them as such, a sequencing artifact in one fragment looks like a
confident variant supported by "many" reads.

The samtools workflow for paired data is a small chain: `collate` (group mate
pairs) → `fixmate -m` (add mate-score tags) → `sort` → `markdup`. Single-end reads
skip collate/fixmate and go straight to `markdup`, which identifies duplicates by
their alignment coordinates.

Crucially, for germline variant calling we **flag** duplicates (leave them in the
BAM with a duplicate bit set) rather than **remove** them — `bcftools mpileup`
ignores flagged duplicates automatically, but keeping them preserves the full record
for inspection. (This differs from the ATAC/ChIP tracks, which remove duplicates
with `markdup -r` because peak-callers prefer a pre-cleaned BAM.)
