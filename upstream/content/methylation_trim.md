# Trimming — fastp (Bisulfite/WGBS)

Bisulfite treatment converts unmethylated cytosines (C) to uracil (read as T), while
methylated cytosines remain C. This means bisulfite reads look unusual — high T content
is expected and normal, not an artifact.

**Extra flags used for bisulfite libraries:**
- `--trim_poly_g` — removes polyG tails common in two-colour Illumina chemistry; these
  appear when sequencing terminates early and the instrument fills in G
- `--length_required 36` — reads shorter than 36 bp after trimming are discarded; Bismark
  needs enough sequence context to distinguish real methylation from conversion artefacts

Standard adapter detection works fine — bisulfite treatment does not affect adapter sequences.
