## Read Trimming — fastp

ChIP-seq reads are trimmed exactly like other assays: fastp removes adapter
sequence and low-quality bases. Clean reads matter here because spurious adapter
or low-quality alignments inflate background and weaken the IP-vs-input contrast
that peak calling depends on.

Paired-end libraries are trimmed as pairs (`--detect_adapter_for_pe`); single-end
libraries use fastp's automatic single-end adapter detection.
