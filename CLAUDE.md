# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

PLACER (Placeability-Aware Long-read Collapsed-Event Reconstruction for Transposable Element insertions) is a bioinformatics tool that detects non-reference TE insertions from long-read BAM alignments (ONT/PacBio). Designed for real WGS BAMs with repetitive genomes, imperfect alignments, and high TE background signal.

## Core Problem

TE insertions in repetitive regions are missed because:
- Mapping is ambiguous (multiple equally good loci)
- Aligners may "linearize" TE sequence (hide split/clip)
- Naive local assembly mixes nearby loci/haplotypes
- Genotyping fails when REF support is ill-defined in repeats

## Architecture (Three-Layer Design)

### Stream Layer (Single-pass BAM streaming)
- Sequential BAM reading, window statistics, read fragment extraction
- Triggers window sealing and enqueues tasks
- Does NOT do POA, EM, or full local alignment

### Task Layer (Async task queue)
- `COMPONENT_BUILD`: bin/subgroup decomposition, candidate locus collection
- `LOCAL_ALIGN`: Restricted local re-alignment
- `ASSEMBLY`: Graph-POA / contig generation
- `COLLAPSE`: Structure-level merging
- `GENOTYPE`: EM (ALT/REF/NULL) + spatial priors

### Output Layer
- `scientific.vcf`: Primary output with locus-set, multi-locus (Tier2), evidence fields
- `engineering.tsv/bed`: Evidence-only output (no representative coordinates by default)

## Key Data Structures

- **ReadSketch**: Lightweight read representation (read_id, chrom, pos, mapq, cigar_ops, SA_summary, probe_fragments)
- **ProbeFragment**: Short TE screening fragments (200-400bp) from END5/END3/SOFTCLIP/CIGAR events/SA breakpoints
- **WindowBuffer**: Sliding window with statistics and buffered reads (triggered when stats exceed quantile thresholds)
- **Component**: Event decomposition unit with locus_set, breakpoint_hypothesis, core_reads
- **StructuralRepresentative**: Structure-level merged output with signature, contig, polymorphism_summary

## Processing Pipeline

1. **Stream traversal**: Single-pass BAM, window statistics, probe extraction, trigger logic
2. **Gate 1**: TE-proxy using k-mer/minimizer matching on probe fragments (cheap screening)
3. **Full TE scan**: Local minimizer chaining on triggered reads only
4. **Component build**: Binning by reference coordinate, locus-set clustering, breakpoint decomposition
5. **Local re-alignment**: Restricted to candidate loci ±10kb (not genome-wide)
6. **Assembly**: Graph-POA (up to 2 output paths per component)
7. **Collapsing**: Structure-level merging by signature
8. **Genotyping**: EM with ALT/REF/NULL mixture model + spatial priors

## Inputs

- `input.bam`: Long-read alignments (ONT/PacBio) with SEQ
- `ref.fa`: Reference genome FASTA (for local re-alignment, not full re-alignment)
- `te.fa`: TE library FASTA (for TE-proxy and annotation)
- Optional: `ref.fa.fai`, `ref.mmi` (for alignment acceleration), MD tag (for mismatch density acceleration)

## Critical Design Constraints

- **Single-pass streaming**: NO full genome scan, NO BAM re-read/seek after initial pass
- **Restricted search spaces**: Local alignment only near candidate loci, not genome-wide
- **Tiered output**:
  - Tier1: Small candidate set, high placeability, single coordinate
  - Tier2: Multiple equivalent loci, consistent structure
  - Tier3: No stable locus or mixed components (evidence-only)

## Implementation Language

C++ (htslib) or Python (pysam) + C++ extensions recommended. Must follow single-pass principle.

## Key Dependencies

- htslib or pysam (BAM handling)
- minimap2 or edlib/ksw2 (local alignment)
- TE library indexing (k-mer/minimizer)

## Implementation Order (MVP to Complete)

1. Stream + WindowBuffer + Trigger (CIGAR/SA/clip only, no TE-proxy yet)
2. Gate 1 (endpoints + clip + I/D probes) + TE k-mer index
3. Component build + restricted candidate loci
4. Local re-alignment (candidate loci ±10kb)
5. POA assembly + structure-level merging
6. Placeability scoring + Tier output
7. EM genotyping + spatial priors
8. TE reverse index (for召回 with missing secondary alignments)
9. Optional: MD peak accelerator

## Core Algorithms

- **TE-proxy**: k-mer/minimizer hit density on probe fragments (not full-read)
- **PlaceabilityScore**: Δ (best vs second-best locus score), SideConsistency, SupportConsistency, CandidateSetSizePenalty
- **EM genotyping**: ALT/REF/NULL mixture with spatial priors, Beta-Binomial posterior intervals
