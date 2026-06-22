## Peak Calling — MACS2 (with input control)

MACS2 identifies genomic regions where ChIP reads pile up above background —
the binding sites of the protein or histone mark that was immunoprecipitated.

**Input control (`-c`).** ChIP-seq's defining feature is the matched **input**: a
sample that skips the antibody pull-down, capturing the background read
distribution (open chromatin, mappability, copy number). MACS2 uses it to decide
which pileups are real enrichment versus background. If a `control` column names an
input sample, it is passed as `-c`; without one, MACS2 models background from the
ChIP sample alone (weaker, but valid for a teaching run).

**Narrow vs broad.** The right peak shape depends on the mark:

- **Narrow** (default): sharp, localized peaks — transcription factors and
  punctate marks like H3K4me3 and H3K27ac. MACS2 builds a fragment-size model
  from the data.
- **Broad** (`--peak-type broad`): wide enrichment domains — repressive or
  gene-body marks like H3K27me3, H3K9me3, H3K36me3. Adds `--broad`, producing a
  `.broadPeak` file.

Note these are **not** the ATAC-seq settings (`--nomodel --shift --extsize`): those
correct the Tn5 transposase cut-site offset, which does not apply to ChIP.
Paired-end libraries are called as fragments (`-f BAMPE`); single-end as `-f BAM`.
