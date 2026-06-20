# Project Brief: HTS Preprocessing Teaching Tool

## For: Claude Code
## Prepared by: Gustavo (instructor, breast cancer immunology lab) for an undergrad course

---

## 1. Purpose

Build a terminal-based, guided teaching tool that lets undergraduate students preprocess
real high-throughput sequencing (HTS) data, step by step, on a shared lab server. Students
should run real bioinformatics tools (not simulations) and get explanations, checkpoints,
and validation feedback as they go. The tool covers three tracks — RNA-seq, ATAC-seq, and
DNA methylation (WGBS/RRBS) — each ending at the point where the instructor's existing
downstream analysis (e.g. DEG analysis for RNA-seq) takes over.

This is NOT a notebook or web app. It is a command-line tool, because:
- It runs on a shared server in a computer lab (not student laptops)
- Command-line fluency is itself part of what students need to learn
- It needs to invoke real CLI bioinformatics tools as subprocesses

---

## 2. Environment & Constraints

- **Host**: shared server in a computer lab, controlled by the instructor (not student
  laptops, not cloud/Colab)
- **Users**: multiple undergrad students with individual accounts/working directories on
  the same server
- **Compute**: assume modest, shared resources — pipelines should run on subsampled
  FASTQ files (small enough to finish within a ~1-3 hour lab session), not full-depth data
- **No GPU assumed** unless confirmed otherwise

---

## 3. Architecture Overview

Three components:

1. **Conda environment** — one shared, pinned environment with all required tools,
   installed once on the server, activated by all students.
2. **CLI teaching tool** (Python package, `typer` + `rich` recommended) — the guided
   walkthrough itself. Installed into the conda environment so it's available as a command
   (e.g. `htsprep`).
3. **Shared data + per-student workspace structure** on the server filesystem.

```
/opt/htsprep/                     # tool install + environment (admin-managed)
/data/htsprep/shared/             # read-only sample FASTQs, reference genomes/indices
/data/htsprep/students/<username>/ # per-student working directories, auto-created
```

---

## 4. Repository Structure (suggested)

```
htsprep/
├── pyproject.toml
├── environment.yml                 # conda spec, pinned versions
├── README.md                       # setup instructions for the instructor/admin
├── htsprep/
│   ├── __init__.py
│   ├── cli.py                      # typer app, entry point
│   ├── tracks/
│   │   ├── shared_qc.py            # Module 0
│   │   ├── rnaseq.py               # Module 1
│   │   ├── atacseq.py              # Module 2
│   │   └── methylation.py          # Module 3
│   ├── checkpoints.py              # validation logic for each step's output
│   ├── content/                    # markdown/text explanation snippets per step
│   │   ├── shared_qc/
│   │   ├── rnaseq/
│   │   ├── atacseq/
│   │   └── methylation/
│   └── progress.py                 # per-student progress logging (for instructor view)
├── scripts/
│   └── prepare_sample_data.sh      # subsamples/downloads reference test datasets
└── tests/
    └── ...
```

---

## 5. Conda Environment — Tools Needed

**Shared / Module 0 (QC fundamentals, all tracks):**
- `fastqc`
- `multiqc`
- `fastp` (trimming, also used per-track)

**RNA-seq track:**
- `star` or `hisat2` (aligner — pick one; STAR is more standard for RNA-seq teaching)
- `samtools`
- `subread` (for `featureCounts`) or `salmon` (alignment-free quantification — simpler,
  may be preferable for undergrads)

**ATAC-seq track:**
- `bowtie2`
- `samtools`
- `picard` (duplicate marking) or `samtools markdup`
- `macs2` (peak calling)
- `deeptools` (optional, for fragment size / signal QC plots)

**Methylation track (WGBS/RRBS):**
- `bismark` (bisulfite-aware alignment + methylation extraction)
- `bowtie2` (Bismark dependency)
- `samtools`
- `trim_galore` (commonly paired with Bismark for bisulfite-aware trimming) — alternative
  to `fastp` for this track specifically; confirm with instructor which trimmer to
  standardize on across all tracks vs. per-track

**CLI tool dependencies:**
- `typer`, `rich`, `python>=3.10`

Claude Code should pin specific versions in `environment.yml` and verify the environment
actually solves/installs cleanly — these bioinformatics tools can have tricky dependency
chains (especially Bismark + Bowtie2 version compatibility).

---

## 6. CLI Tool Behavior Spec

### Command structure (proposed, adjust as needed)
```
htsprep list                          # show available tracks/modules
htsprep start <track>                 # initializes student's working dir for a track
htsprep step <track> <step-name>      # prints explanation, runs/guides the step
htsprep check <track> <step-name>     # validates output before unlocking next step
htsprep progress                      # shows student's own progress
htsprep progress --all                # instructor-only: shows all students' progress
```

### Per-step behavior
Each step should:
1. Print a short explanation (markdown rendered via `rich`) of *why* this step matters
   and what's about to happen
2. Either (a) run the real command for the student and show live output, or
   (b) print the exact command and have the student run it themselves, then validate —
   **decision needed**: fully automated vs. "type the command yourself" pedagogy.
   (Recommendation: automate Module 0 to build confidence, then shift to "type it
   yourself with guidance" for later modules so students build real CLI fluency.)
3. Validate output before allowing progression (e.g. checking that a FastQC report was
   generated, an alignment BAM is non-empty and has a sane mapping rate, a peak file has
   a reasonable number of peaks, etc.)
4. Log completion + timestamp to a per-student progress file for instructor visibility
5. On failure, give a specific, actionable error — not just "step failed"

### Checkpoints / understanding checks
At key points (e.g. after QC, after alignment), prompt the student with a short question
(e.g. "What does a mapping rate below 70% suggest?") before unlocking the next step.
**Decision needed**: free-text logged for instructor review, or multiple-choice
auto-graded, or just a "did you discuss this with your TA" gate.

---

## 7. Pipeline Content Per Track

### Module 0 — Shared QC fundamentals (prerequisite for all tracks)
- Raw read inspection (`fastqc`)
- Phred quality scores, what they mean
- Adapter contamination, GC content skew, duplication levels — what's normal vs.
  concerning
- `multiqc` to aggregate reports
- Ends with: students can read a FastQC/MultiQC report and identify a problem

### Module 1 — RNA-seq
- Trimming (`fastp`)
- Alignment (`STAR` or `HISAT2` to a reference genome/transcriptome)
- Quantification (`featureCounts` or `salmon`) → gene-level counts matrix
- **Ends at**: a counts matrix ready for the instructor's existing DESeq2/edgeR workflow

### Module 2 — ATAC-seq
- Trimming (`fastp`)
- Alignment (`Bowtie2`)
- Filtering: remove mitochondrial reads, mark/remove duplicates, fragment size QC
  (nucleosome-free vs. mono-nucleosome fragments)
- Peak calling (`MACS2`)
- **Ends at**: a peak set / accessibility matrix

### Module 3 — Methylation (WGBS/RRBS)
- Bisulfite-aware trimming (`Trim Galore` or `fastp` with appropriate settings)
- Bisulfite alignment (`Bismark`)
- Methylation extraction (`bismark_methylation_extractor`)
- **Ends at**: per-CpG methylation calls, ready for downstream DMR analysis

---

## 8. Sample Data Requirements

Need small, public, subsampled datasets for each track that:
- Are real (not synthetic) so QC issues/patterns are genuine
- Are small enough to process in a single lab session on shared/modest compute
- Have a known "good" outcome so checkpoint validation has a ground truth to check against

**Open task for Claude Code / instructor**: identify specific accessions (e.g. from SRA,
ENCODE, GEO) for each track and write a `scripts/prepare_sample_data.sh` that downloads
and subsamples them (e.g. via `seqtk sample`) to a manageable read count. This brief
intentionally does not hardcode accessions — they should be selected and verified for
appropriateness (organism matching the course's focus, file size, known biology) at
build time.

**Optional, given lab's TME/cytokine focus**: if helpful for student engagement, sample
datasets could be chosen from breast cancer or immune cell contexts (matching the
instructor's own research area), but this isn't a requirement — any well-behaved public
dataset per assay type works.

---

## 9. Per-Student State & Instructor Visibility

- Each student's progress, working files, and logs live under
  `/data/htsprep/students/<username>/`
- A simple progress log (JSON or SQLite — Claude Code's choice) records: track, step,
  timestamp, pass/fail, and any checkpoint answers
- Instructor command (`htsprep progress --all`) aggregates this across students —
  needs a way to distinguish instructor/TA accounts from student accounts (e.g. a config
  file listing admin usernames, or Unix group membership)

---

## 10. Open Decisions for Claude Code to Flag/Resolve With the Instructor

1. Aligner choice for RNA-seq: STAR vs. HISAT2; quantification via featureCounts vs.
   Salmon
2. Trimmer standardization: one trimmer (fastp) across all tracks, or Trim Galore
   specifically for the methylation track
3. Step automation level: fully automated vs. "student types the command" — likely
   should differ by module (see Section 6)
4. Checkpoint question format: free text, multiple choice, or instructor-review gate
5. Specific sample dataset accessions per track (Section 8)
6. Reference genome(s) needed (organism, genome build) — affects index sizes and storage
   on the server
7. Authentication/admin distinction for progress visibility (Unix groups vs. config file)

---

## 11. Final Output Format — Must Match OBAMA Pipeline Input

All tracks' final output (the file handed off to downstream analysis) must conform to the
input format required by the OBAMA pipeline (https://github.com/AOG-Lab/OBAMA), regardless
of which track produced it (RNA-seq counts, ATAC-seq accessibility matrix, methylation
calls).

**Required structure:**
- File format: `.csv` or `.tsv`
- Each **row** = one sample
- Each **column** = one gene (or feature, depending on track)
- The **first two columns must be, in this order**:
  1. `geo_accession` — unique sample identifier
  2. `disease.state` — biological condition label

**Critical naming rule (important even when it seems unintuitive):**
OBAMA's code only reads columns literally named `control` and `disease` for its group
comparisons. **Even if a given dataset doesn't actually have a disease/control
distinction** (e.g. a time-course, a dose series, or some other two-group comparison),
the two groups must still be labeled `control` and `disease` in the `disease.state`
column — otherwise OBAMA will not recognize the grouping at all. Claude Code should bake
this naming requirement into the final-output-writing step for every track, with a clear
comment in the code explaining *why* the labels are forced this way (so future
maintainers don't "fix" it into something more semantically accurate and break
compatibility with OBAMA).

Example of the required final structure:

| geo_accession | disease.state | GENE1 | GENE2 | ... |
|---|---|---|---|---|
| GSM650656 | disease | 192.72 | 94.82 | ... |
| GSM650657 | control | 241.33 | 120.10 | ... |

This should be validated as part of each track's final `check` step — i.e. the checkpoint
for the last step of every module should confirm the output file has exactly this
structure before marking the module complete.

---

## 12. Success Criteria

- A student with no prior CLI bioinformatics experience can complete one full track in
  a single lab session
- Every step either runs a real tool or has the student run a real tool — no simulated
  output
- The instructor can see, after the fact, where each student is/was stuck
- The RNA-seq track's final output drops directly into the instructor's existing
  DESeq2/edgeR workflow with no reformatting needed
