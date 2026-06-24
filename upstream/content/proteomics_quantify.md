# Proteomics — From Spectra to a Protein Matrix

**The full workflow.** In bottom-up proteomics, proteins are digested into peptides,
separated by liquid chromatography, and measured by a mass spectrometer (LC-MS/MS).
A search engine — most commonly **MaxQuant** (with its Andromeda search) or
**FragPipe** (MSFragger) — matches the spectra to a protein database, controls the
false-discovery rate, and quantifies each protein. Its main output table is
**`proteinGroups.txt`**.

> **This step is upstream and out of scope.** Running MaxQuant/FragPipe on raw
> spectra (`.raw`/`.mzML`) needs the instrument files and heavy, mostly-GUI tools.
> This track starts from the **intensity matrix** those tools produce — exactly as
> the methylation array track starts from a GEO beta matrix rather than raw IDATs.
> The analysis below is real; the search that produced your input is assumed done.

**What this track does** (real, reproducible computation):

1. **Drop non-protein rows** — entries flagged as contaminants, decoy (`Reverse`)
   hits, or "only identified by site" are removed.
2. **Filter by valid values** — a protein is kept only if it was actually quantified
   in enough samples *within each group* (the `--min-valid` fraction, default 0.5).
   Missing values in proteomics are informative — a protein may be genuinely absent
   (below the detection limit), not just unmeasured — so we filter rather than trust
   sparse rows.
3. **log2-transform** — intensities span several orders of magnitude; log2 makes the
   distribution roughly symmetric and the differences additive (a log2 difference of
   1 = a 2-fold change), which is what limma models.
4. **Normalize** — `median` centering (default) aligns each sample's median intensity
   to correct for differences in total loaded protein, without forcing the
   distributions to be identical. `quantile` is stricter (makes all samples share one
   distribution); `none` defers normalization to R.

**Quantification flavors.** **LFQ** (label-free quantification, the default) is
normalized for between-sample comparison; **iBAQ** is roughly proportional to molar
abundance within a sample; raw **Intensity** is unnormalized. Pick with
`--intensity-col`.

We do **not** impute the remaining missing values — imputation (e.g.
downshifted-normal) is a modeling choice best made deliberately in your analysis, and
limma handles missing values.
