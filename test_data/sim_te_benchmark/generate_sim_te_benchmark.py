#!/usr/bin/env python3
"""Generate synthetic ONT TE benchmark with insertion sizes in [300, 2000] bp.

TE templates are real curated consensus sequences from Dfam:
- L1HS_5end (DF000000226, v4)
- L1HS_3end (DF000000225, v4)
- AluYa5 (DF000000053, v4)
- SVA_F (DF000001072, v4)
"""

from __future__ import annotations

import argparse
import random
import subprocess
from pathlib import Path
from typing import Dict, List, Tuple

CURATED_TEMPLATE_BASENAME = "te_consensus_dfam.fa"


def random_dna(rng: random.Random, length: int) -> str:
    return "".join(rng.choice("ACGT") for _ in range(length))


def read_single_fasta(path: Path) -> Tuple[str, str]:
    header = ""
    seq_chunks: List[str] = []
    with path.open() as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header:
                    raise ValueError(f"{path} has more than one FASTA record")
                header = line[1:]
            else:
                seq_chunks.append(line.upper())
    if not header:
        raise ValueError(f"{path} has no FASTA record")
    return header, "".join(seq_chunks)


def read_multi_fasta(path: Path) -> Dict[str, str]:
    records: Dict[str, str] = {}
    header = ""
    seq_chunks: List[str] = []

    def flush() -> None:
        nonlocal header, seq_chunks
        if not header:
            return
        token = header.split()[0]
        seq = "".join(seq_chunks).upper()
        if token in records:
            raise ValueError(f"duplicate FASTA record: {token}")
        records[token] = seq
        header = ""
        seq_chunks = []

    with path.open() as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                flush()
                header = line[1:]
            else:
                seq_chunks.append(line)
    flush()
    if not records:
        raise ValueError(f"{path} has no FASTA records")
    return records


def sanitize_consensus(seq: str) -> str:
    # Dfam consensus can contain ambiguity codes; replace with A for deterministic A/C/G/T-only templates.
    return "".join(c if c in "ACGT" else "A" for c in seq.upper())


def build_curated_te_library(consensus_fasta: Path) -> Dict[str, str]:
    recs = read_multi_fasta(consensus_fasta)
    required = {
        "L1HS_5end_DF000000226_v4",
        "L1HS_3end_DF000000225_v4",
        "AluYa5_DF000000053_v4",
        "SVA_F_DF000001072_v4",
    }
    missing = sorted(required.difference(recs.keys()))
    if missing:
        raise ValueError(
            f"missing required consensus records in {consensus_fasta}: {', '.join(missing)}"
        )

    # L1HS in Dfam is split as 5' and 3' models. Concatenate both curated parts to build a long template.
    l1 = sanitize_consensus(recs["L1HS_5end_DF000000226_v4"] + recs["L1HS_3end_DF000000225_v4"])
    alu = sanitize_consensus(recs["AluYa5_DF000000053_v4"])
    sva = sanitize_consensus(recs["SVA_F_DF000001072_v4"])

    return {
        "L1:L1HS": l1,
        "Alu:AluYa5": alu,
        "SVA:SVA_F": sva,
    }


def write_fasta(path: Path, records: List[Tuple[str, str]], line_width: int = 80) -> None:
    with path.open("w") as out:
        for name, seq in records:
            out.write(f">{name}\n")
            for i in range(0, len(seq), line_width):
                out.write(seq[i : i + line_width] + "\n")


def clamp_start(start: int, min_start: int, max_start: int) -> int:
    return max(min_start, min(start, max_start))


def sam_record(
    qname: str,
    flag: int,
    chrom: str,
    pos: int,
    mapq: int,
    cigar: str,
    seq: str,
    tags: List[str] | None = None,
) -> str:
    qual = "I" * len(seq)
    fields = [
        qname,
        str(flag),
        chrom,
        str(pos),
        str(mapq),
        cigar,
        "*",
        "0",
        "0",
        seq,
        qual,
    ]
    if tags:
        fields.extend(tags)
    return "\t".join(fields)


def write_tsv(path: Path, header: List[str], rows: List[List[str]]) -> None:
    with path.open("w") as out:
        out.write("\t".join(header) + "\n")
        for row in rows:
            out.write("\t".join(map(str, row)) + "\n")


def build_event_table() -> List[Dict[str, object]]:
    return [
        {
            "event_id": "TP_L1_INS_STRONG",
            "chrom": "chrTEST",
            "pos": 32000,
            "family": "L1",
            "truth_label": "POS",
            "challenge_type": "BASELINE",
            "expected_call": "YES",
            "scenario": "long_insertion_I_2000bp",
            "support_reads": 12,
            "note": "Canonical long L1 insertion support",
            "kind": "ins",
            "ins_source": "L1",
            "ins_len": 2000,
            "polyA_len": 8,
            "mapq": 60,
        },
        {
            "event_id": "TP_ALU_INS_STRONG",
            "chrom": "chrTEST",
            "pos": 76000,
            "family": "Alu",
            "truth_label": "POS",
            "challenge_type": "BASELINE",
            "expected_call": "YES",
            "scenario": "long_insertion_I_310bp",
            "support_reads": 10,
            "note": "Canonical AluYa5 insertion support",
            "kind": "ins",
            "ins_source": "Alu",
            "ins_len": 310,
            "polyA_len": 20,
            "mapq": 60,
        },
        {
            "event_id": "TP_L1_SPLIT_CLIP",
            "chrom": "chrTEST",
            "pos": 108000,
            "family": "L1",
            "truth_label": "POS",
            "challenge_type": "BASELINE",
            "expected_call": "YES",
            "scenario": "split_plus_long_softclip_500bp",
            "support_reads": 6,
            "note": "Primary + supplementary split-read support",
            "kind": "split",
            "split_source": "L1",
            "clip_len": 500,
            "clip_polyA_len": 25,
            "match_len": 600,
            "mapq": 60,
        },
        {
            "event_id": "F1_FP_RANDOM_INS",
            "chrom": "chrTEST",
            "pos": 138000,
            "family": "NA",
            "truth_label": "NEG",
            "challenge_type": "F1_FP",
            "expected_call": "NO",
            "scenario": "non_te_long_insertion_900bp",
            "support_reads": 10,
            "note": "Large insertion but sequence is non-TE random DNA",
            "kind": "ins",
            "ins_source": "RANDOM",
            "ins_len": 900,
            "mapq": 60,
        },
        {
            "event_id": "F1_FP_POLYA_SOFTCLIP",
            "chrom": "chrTEST",
            "pos": 166000,
            "family": "NA",
            "truth_label": "NEG",
            "challenge_type": "F1_FP",
            "expected_call": "NO",
            "scenario": "polyA_softclip_only_350bp",
            "support_reads": 8,
            "note": "PolyA-rich soft-clip that can mimic TE tails",
            "kind": "softclip",
            "clip_len": 350,
            "match_len": 600,
            "mapq": 60,
        },
        {
            "event_id": "F2_FN_LOW_SUPPORT_L1",
            "chrom": "chrTEST",
            "pos": 184000,
            "family": "L1",
            "truth_label": "POS",
            "challenge_type": "F2_FN",
            "expected_call": "NO",
            "scenario": "te_insertion_low_support_1600bp",
            "support_reads": 2,
            "note": "Real TE insertion but support below target threshold",
            "kind": "ins",
            "ins_source": "L1",
            "ins_len": 1600,
            "polyA_len": 8,
            "mapq": 60,
        },
        {
            "event_id": "F2_FN_SHORT_INS_ALU",
            "chrom": "chrTEST",
            "pos": 193000,
            "family": "Alu",
            "truth_label": "POS",
            "challenge_type": "F2_FN",
            "expected_call": "NO",
            "scenario": "te_insertion_low_end_300bp",
            "support_reads": 7,
            "note": "Lower-end long insertion (300bp) to test boundary behavior",
            "kind": "ins",
            "ins_source": "Alu",
            "ins_len": 300,
            "polyA_len": 20,
            "mapq": 60,
        },
        {
            "event_id": "TP_L1_INS_SECONDARY",
            "chrom": "chrTEST",
            "pos": 68000,
            "family": "L1",
            "truth_label": "POS",
            "challenge_type": "BASELINE",
            "expected_call": "YES",
            "scenario": "long_insertion_I_2000bp_locus2",
            "support_reads": 12,
            "note": "Second long L1 insertion in an alternate locus",
            "kind": "ins",
            "ins_source": "L1",
            "ins_len": 2000,
            "polyA_len": 8,
            "mapq": 60,
        },
    ]


def take_prefix_with_a_pad(seq: str, length: int) -> str:
    if length <= len(seq):
        return seq[:length]
    return seq + ("A" * (length - len(seq)))


def take_suffix_with_a_pad(seq: str, length: int) -> str:
    if length <= len(seq):
        return seq[-length:]
    return ("A" * (length - len(seq))) + seq


def build_te_like_insert(source: str, ins_len: int, poly_a_len: int, te_lib: Dict[str, str]) -> str:
    if ins_len <= 0:
        raise ValueError("ins_len must be positive")
    poly_a_len = max(0, min(poly_a_len, max(0, ins_len - 40)))
    core_len = ins_len - poly_a_len

    if source == "L1":
        # Use a higher-complexity L1 region to keep long ONT remapping stable in TLDR.
        return take_prefix_with_a_pad(te_lib["L1:L1HS"], core_len) + ("A" * poly_a_len)
    if source == "SVA":
        # Use a higher-complexity SVA region to keep long ONT remapping stable in TLDR.
        return take_prefix_with_a_pad(te_lib["SVA:SVA_F"], core_len) + ("A" * poly_a_len)
    if source == "Alu":
        # Alu insertions are near full-length; use consensus prefix with tail extension if needed.
        return take_prefix_with_a_pad(te_lib["Alu:AluYa5"], core_len) + ("A" * poly_a_len)
    raise ValueError(f"unknown TE insertion source: {source}")


def build_insert_seq(event: Dict[str, object], te_lib: Dict[str, str], rng: random.Random) -> str:
    source = str(event.get("ins_source"))
    ins_len = int(event["ins_len"])
    poly_a_len = int(event.get("polyA_len", 0))
    if source == "RANDOM":
        return random_dna(rng, ins_len)
    return build_te_like_insert(source, ins_len, poly_a_len, te_lib)


def main() -> None:
    parser = argparse.ArgumentParser(description="Generate sim_te_benchmark test set")
    parser.add_argument(
        "--outdir",
        default=str(Path(__file__).resolve().parent),
        help="output directory (default: test_data/sim_te_benchmark)",
    )
    parser.add_argument(
        "--te-consensus-fasta",
        default=str(Path(__file__).resolve().parent / CURATED_TEMPLATE_BASENAME),
        help="curated TE consensus FASTA (default: test_data/sim_te_benchmark/te_consensus_dfam.fa)",
    )
    parser.add_argument("--seed", type=int, default=20260218, help="random seed")
    args = parser.parse_args()

    outdir = Path(args.outdir).resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    consensus_fasta = Path(args.te_consensus_fasta).resolve()

    ref_name, ref_seq = read_single_fasta(outdir / "ref.fa")
    ref_len = len(ref_seq)

    rng = random.Random(args.seed)
    events = build_event_table()

    te_lib = build_curated_te_library(consensus_fasta)
    write_fasta(outdir / "te_tldr.fa", list(te_lib.items()))

    manifest_rows: List[List[str]] = []
    sam_lines: List[str] = []

    # Background reads (long ONT-like aligned reads).
    background_count = 140
    background_len = 400
    bg_start_max = ref_len - background_len + 1
    for idx in range(background_count):
        qname = f"BG_{idx:04d}"
        pos = rng.randint(1, bg_start_max)
        seq = ref_seq[pos - 1 : pos - 1 + background_len]
        cigar = f"{background_len}M"
        sam_lines.append(sam_record(qname, 0, ref_name, pos, 60, cigar, seq))
        manifest_rows.append([qname, "0", ref_name, str(pos), "60", cigar])

    left_flank = 600
    right_flank = 600

    for event in events:
        event_id = str(event["event_id"])
        chrom = str(event["chrom"])
        pos = int(event["pos"])
        support_reads = int(event["support_reads"])
        mapq = int(event["mapq"])
        kind = str(event["kind"])

        if kind == "ins":
            insert_seq = build_insert_seq(event, te_lib, rng)
            ins_len = len(insert_seq)
            for idx in range(support_reads):
                jitter = rng.randint(-3, 3)
                start = pos - left_flank + jitter
                start = clamp_start(start, 1, ref_len - left_flank - right_flank + 1)
                left_ref = ref_seq[start - 1 : start - 1 + left_flank]
                right_ref = ref_seq[start - 1 + left_flank : start - 1 + left_flank + right_flank]
                read_seq = left_ref + insert_seq + right_ref
                cigar = f"{left_flank}M{ins_len}I{right_flank}M"
                qname = f"{event_id}_INS_{idx:03d}"
                sam_lines.append(sam_record(qname, 0, chrom, start, mapq, cigar, read_seq))
                manifest_rows.append([qname, "0", chrom, str(start), str(mapq), cigar])

        elif kind == "softclip":
            clip_len = int(event["clip_len"])
            match_len = int(event["match_len"])
            for idx in range(support_reads):
                jitter = rng.randint(-3, 3)
                start = pos - match_len + jitter
                start = clamp_start(start, 1, ref_len - match_len + 1)
                left_ref = ref_seq[start - 1 : start - 1 + match_len]
                read_seq = left_ref + ("A" * clip_len)
                cigar = f"{match_len}M{clip_len}S"
                qname = f"{event_id}_SC3_{idx:03d}"
                sam_lines.append(sam_record(qname, 0, chrom, start, mapq, cigar, read_seq))
                manifest_rows.append([qname, "0", chrom, str(start), str(mapq), cigar])

        elif kind == "split":
            clip_len = int(event["clip_len"])
            match_len = int(event["match_len"])
            clip_poly_a_len = int(event.get("clip_polyA_len", 0))
            clip_seq = build_te_like_insert("L1", clip_len, clip_poly_a_len, te_lib)
            for idx in range(support_reads):
                p_jitter = rng.randint(-3, 3)
                s_jitter = rng.randint(-3, 3)
                primary_start = pos - match_len + p_jitter
                primary_start = clamp_start(primary_start, 1, ref_len - match_len + 1)
                supp_start = pos + 40 + s_jitter
                supp_start = clamp_start(supp_start, 1, ref_len - match_len + 1)

                p_ref = ref_seq[primary_start - 1 : primary_start - 1 + match_len]
                s_ref = ref_seq[supp_start - 1 : supp_start - 1 + match_len]
                p_seq = p_ref + clip_seq
                s_seq = clip_seq + s_ref
                p_cigar = f"{match_len}M{clip_len}S"
                s_cigar = f"{clip_len}S{match_len}M"
                qname = f"{event_id}_SPLIT_{idx:03d}"

                p_sa = f"SA:Z:{chrom},{supp_start},+,{s_cigar},{mapq},0;"
                s_sa = f"SA:Z:{chrom},{primary_start},+,{p_cigar},{mapq},0;"
                sam_lines.append(sam_record(qname, 0, chrom, primary_start, mapq, p_cigar, p_seq, [p_sa]))
                sam_lines.append(sam_record(qname, 2048, chrom, supp_start, mapq, s_cigar, s_seq, [s_sa]))
                manifest_rows.append([qname, "0", chrom, str(primary_start), str(mapq), p_cigar])
                manifest_rows.append([qname, "2048", chrom, str(supp_start), str(mapq), s_cigar])

        else:
            raise ValueError(f"unknown event kind: {kind}")

    # Write read manifest.
    write_tsv(
        outdir / "read_manifest.tsv",
        ["qname", "flag", "chrom", "pos", "mapq", "cigar"],
        manifest_rows,
    )

    # Write truth tables.
    truth_rows = [
        [
            str(e["event_id"]),
            str(e["chrom"]),
            str(e["pos"]),
            str(e["family"]),
            str(e["truth_label"]),
            str(e["challenge_type"]),
            str(e["expected_call"]),
            str(e["scenario"]),
            str(e["support_reads"]),
            str(e["note"]),
        ]
        for e in events
    ]
    write_tsv(
        outdir / "truth_events.tsv",
        [
            "event_id",
            "chrom",
            "pos",
            "family",
            "truth_label",
            "challenge_type",
            "expected_call",
            "scenario",
            "support_reads",
            "note",
        ],
        truth_rows,
    )

    truth_pos_strong_ids = {"TP_L1_INS_STRONG", "TP_ALU_INS_STRONG", "TP_L1_SPLIT_CLIP"}
    truth_pos_all_rows = []
    truth_pos_strong_rows = []
    truth_neg_rows = []
    eval_rows = []
    for e in events:
        row = [str(e["chrom"]), str(e["pos"]), str(e["family"]), str(e["event_id"])]
        if e["truth_label"] == "POS":
            truth_pos_all_rows.append(row)
            if e["event_id"] in truth_pos_strong_ids:
                truth_pos_strong_rows.append(row)
        else:
            truth_neg_rows.append(row)

        eval_rows.append(
            [
                str(e["event_id"]),
                str(e["chrom"]),
                str(e["pos"]),
                str(e["family"]),
                str(e["truth_label"]),
                str(e["challenge_type"]),
                str(e["expected_call"]),
                "0",
                "NA",
            ]
        )

    write_tsv(
        outdir / "truth_positive_strong.tsv",
        ["chrom", "pos", "family", "event_id"],
        truth_pos_strong_rows,
    )
    write_tsv(
        outdir / "truth_positive_all.tsv",
        ["chrom", "pos", "family", "event_id"],
        truth_pos_all_rows,
    )
    write_tsv(
        outdir / "truth_negative_fp.tsv",
        ["chrom", "pos", "family", "event_id"],
        truth_neg_rows,
    )
    write_tsv(
        outdir / "eval_matches.tsv",
        [
            "event_id",
            "chrom",
            "pos",
            "family",
            "truth_label",
            "challenge_type",
            "expected_call",
            "matched_calls",
            "matched_call_ids",
        ],
        eval_rows,
    )

    # Write SAM.
    sam_path = outdir / "sim_te_benchmark.sam"
    with sam_path.open("w") as sam:
        sam.write("@HD\tVN:1.6\tSO:unsorted\n")
        sam.write(f"@SQ\tSN:{ref_name}\tLN:{ref_len}\n")
        for line in sam_lines:
            sam.write(line + "\n")

    unsorted_bam = outdir / "sim_te_benchmark.unsorted.bam"
    sorted_bam = outdir / "sim_te_benchmark.bam"
    subprocess.run(
        ["samtools", "view", "-b", "-o", str(unsorted_bam), str(sam_path)],
        check=True,
    )
    subprocess.run(
        ["samtools", "sort", "-o", str(sorted_bam), str(unsorted_bam)],
        check=True,
    )
    subprocess.run(["samtools", "index", str(sorted_bam)], check=True)

    # Keep only final BAM + BAI + manifest/truth tables.
    if unsorted_bam.exists():
        unsorted_bam.unlink()
    if sam_path.exists():
        sam_path.unlink()


if __name__ == "__main__":
    main()
