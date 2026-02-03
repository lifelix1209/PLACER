# What is PLACER

PLACER detects non-reference transposable element (TE) insertions from long-read BAM alignments (ONT/PacBio). It is designed for real WGS BAMs: repetitive genomes, imperfect alignments, missing secondary records, and high background TE signal.

PLACER stands for:

Placeability-Aware Long-read Collapsed-Event Reconstruction for Transposable Element insertions (BAM input)

# What problem it solves

TE insertion scalling breaks in repetitive regions because:

mapping is ambiguous (multiple equally good loci)

aligners may “linearize” TE sequence (hide split/clip)

naive local assembly mixes nearby loci / haplotypes

genotyping fails when REF support is ill-defined in repeats

# PLACER addresses this by:

explicitly modeling placeability (unique vs ambiguous placement)

using single-pass streaming over BAM (I/O is the bottleneck)

extracting short probe fragments for cheap TE screening (Gate 1)

reconstructing candidate insertion sequences only where warranted

genotyping with an ALT/REF/NULL mixture model plus spatial priors