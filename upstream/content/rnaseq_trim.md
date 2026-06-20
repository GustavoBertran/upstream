# Adapter Trimming — fastp (RNA-seq)

Two problems hide in raw RNA-seq reads:

**Adapter contamination** — if an RNA fragment is shorter than ~150 bp, the sequencer reads
through the insert and into the synthetic adapter. Adapters do not exist in the genome, so
reads containing them will fail to align or map to wrong locations.

**Low-quality 3′ ends** — Phred quality scores drop toward the end of reads. Bases with
Q < 20 have a >1% error rate; keeping them inflates false-positive variant calls.

`fastp` solves both in a single pass and writes a JSON report used to confirm the step worked.
Look for > 90% of reads passing filters in the fastp output.
