#!/usr/bin/env python3
"""BLAST an insertion sequence against a TE library and report similarity."""

from __future__ import annotations

import argparse
import csv
import hashlib
import shutil
import subprocess
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import List, NoReturn, Sequence, TextIO


REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_TE_LIBRARY = (
    REPO_ROOT
    / "test_data"
    / "MWCichlidTE-3.2.splitted.namefixed_renamed_all_all.fa"
)
IUPAC_DNA = set("ACGTRYSWKMBDHVN")


@dataclass(frozen=True)
class BlastHit:
    query_id: str
    subject_id: str
    percent_identity: float
    alignment_length: int
    mismatches: int
    gap_opens: int
    query_start: int
    query_end: int
    subject_start: int
    subject_end: int
    evalue: str
    bit_score: float
    query_length: int
    subject_length: int
    query_coverage_pct: float


def fail(message: str) -> NoReturn:
    raise SystemExit(message)


def normalize_sequence(text: str) -> str:
    seq = []
    for raw_ch in text:
        ch = raw_ch.upper()
        if ch in {" ", "\t", "\r", "\n", "-"}:
            continue
        if ch == "U":
            ch = "T"
        if ch not in IUPAC_DNA:
            fail(f"invalid sequence character: {raw_ch!r}")
        seq.append(ch)
    normalized = "".join(seq)
    if not normalized:
        fail("empty insertion sequence")
    return normalized


def sequence_from_text(text: str) -> str:
    lines = []
    for line in text.splitlines():
        if line.startswith(">"):
            continue
        lines.append(line)
    return normalize_sequence("".join(lines))


def load_sequence_from_file(path: Path) -> str:
    if not path.is_file():
        fail(f"sequence file not found: {path}")
    return sequence_from_text(path.read_text(encoding="utf-8"))


def wrap_fasta_sequence(seq: str, width: int = 80) -> str:
    chunks = [seq[i : i + width] for i in range(0, len(seq), width)]
    return ">insertion\n" + "\n".join(chunks) + "\n"


def parse_args(argv: Sequence[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Compare an insertion sequence with the TE library using BLASTN and "
            "report identity/query-coverage similarity."
        )
    )
    input_group = parser.add_mutually_exclusive_group()
    input_group.add_argument(
        "--sequence",
        help="Insertion sequence as raw bases. FASTA headers are also accepted.",
    )
    input_group.add_argument(
        "--sequence-file",
        type=Path,
        help="File containing the insertion sequence, either raw bases or FASTA.",
    )
    parser.add_argument(
        "--te-library",
        type=Path,
        default=DEFAULT_TE_LIBRARY,
        help=f"TE library FASTA. Default: {DEFAULT_TE_LIBRARY}",
    )
    parser.add_argument(
        "--db-dir",
        type=Path,
        default=Path("placer_out/te_blastdb"),
        help="Directory for the generated BLAST database.",
    )
    parser.add_argument(
        "--top",
        type=int,
        default=10,
        help="Number of BLAST hits to report.",
    )
    parser.add_argument(
        "--task",
        default="blastn",
        choices=["blastn", "blastn-short", "megablast", "dc-megablast"],
        help="BLASTN task. Use blastn-short for very short insertions.",
    )
    parser.add_argument(
        "--evalue",
        default="10",
        help="BLAST e-value threshold passed to blastn.",
    )
    parser.add_argument(
        "--blastn",
        default="blastn",
        help="blastn executable path or name.",
    )
    parser.add_argument(
        "--makeblastdb",
        default="makeblastdb",
        help="makeblastdb executable path or name.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help="Optional TSV output path. Default: stdout.",
    )
    return parser.parse_args(argv)


def require_file(path: Path, label: str) -> Path:
    resolved = path.expanduser().resolve()
    if not resolved.is_file():
        fail(f"{label} not found: {resolved}")
    return resolved


def require_executable(name_or_path: str, label: str) -> str:
    path = Path(name_or_path).expanduser()
    if path.parent != Path("."):
        resolved = path.resolve()
        if resolved.is_file() and resolved.stat().st_mode & 0o111:
            return str(resolved)
        fail(f"{label} is not executable: {resolved}")
    found = shutil.which(name_or_path)
    if not found:
        local_tool = REPO_ROOT / ".tools" / "blastplus" / "bin" / name_or_path
        if local_tool.is_file() and local_tool.stat().st_mode & 0o111:
            return str(local_tool)
        fail(
            f"{label} not found in PATH: {name_or_path}. "
            "Install NCBI BLAST+ or pass the executable path explicitly."
        )
    return found


def te_library_digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()[:12]


def blast_db_exists(prefix: Path) -> bool:
    return any(prefix.parent.glob(prefix.name + ".n*"))


def build_blast_db(
    *,
    te_library: Path,
    db_dir: Path,
    makeblastdb_exe: str,
) -> Path:
    digest = te_library_digest(te_library)
    db_dir = db_dir.expanduser().resolve() / digest
    db_dir.mkdir(parents=True, exist_ok=True)
    prefix = db_dir / "te_library"
    if blast_db_exists(prefix):
        return prefix

    completed = subprocess.run(
        [
            makeblastdb_exe,
            "-in",
            str(te_library),
            "-dbtype",
            "nucl",
            "-out",
            str(prefix),
        ],
        check=False,
        text=True,
        capture_output=True,
    )
    if completed.returncode != 0:
        fail(
            "makeblastdb failed\n"
            + completed.stdout
            + ("\n" if completed.stdout and completed.stderr else "")
            + completed.stderr
        )
    return prefix


def parse_blast_tabular(text: str) -> List[BlastHit]:
    hits: List[BlastHit] = []
    for line in text.splitlines():
        if not line.strip():
            continue
        fields = line.rstrip("\n").split("\t")
        if len(fields) != 14:
            fail(f"malformed BLAST tabular row: {line}")
        qseqid, sseqid = fields[0], fields[1]
        pident = float(fields[2])
        length = int(fields[3])
        mismatch = int(fields[4])
        gapopen = int(fields[5])
        qstart = int(fields[6])
        qend = int(fields[7])
        sstart = int(fields[8])
        send = int(fields[9])
        evalue = fields[10]
        bitscore = float(fields[11])
        qlen = int(fields[12])
        slen = int(fields[13])
        query_span = abs(qend - qstart) + 1
        qcov = (query_span / qlen * 100.0) if qlen > 0 else 0.0
        hits.append(
            BlastHit(
                query_id=qseqid,
                subject_id=sseqid,
                percent_identity=pident,
                alignment_length=length,
                mismatches=mismatch,
                gap_opens=gapopen,
                query_start=qstart,
                query_end=qend,
                subject_start=sstart,
                subject_end=send,
                evalue=evalue,
                bit_score=bitscore,
                query_length=qlen,
                subject_length=slen,
                query_coverage_pct=qcov,
            )
        )
    return hits


def run_blast(
    *,
    seq: str,
    db_prefix: Path,
    blastn_exe: str,
    task: str,
    evalue: str,
    top: int,
) -> List[BlastHit]:
    if top <= 0:
        fail("--top must be > 0")
    with tempfile.TemporaryDirectory() as tmpdir:
        query_path = Path(tmpdir) / "insertion.fa"
        query_path.write_text(wrap_fasta_sequence(seq), encoding="utf-8")
        completed = subprocess.run(
            [
                blastn_exe,
                "-task",
                task,
                "-query",
                str(query_path),
                "-db",
                str(db_prefix),
                "-evalue",
                evalue,
                "-max_target_seqs",
                str(top),
                "-outfmt",
                "6 qseqid sseqid pident length mismatch gapopen qstart qend "
                "sstart send evalue bitscore qlen slen",
            ],
            check=False,
            text=True,
            capture_output=True,
        )
    if completed.returncode != 0:
        fail(
            "blastn failed\n"
            + completed.stdout
            + ("\n" if completed.stdout and completed.stderr else "")
            + completed.stderr
        )
    return parse_blast_tabular(completed.stdout)


def write_hits_tsv(handle: TextIO, hits: Sequence[BlastHit]) -> None:
    fieldnames = [
        "rank",
        "te_id",
        "percent_identity",
        "query_coverage_pct",
        "alignment_length",
        "query_length",
        "subject_length",
        "mismatches",
        "gap_opens",
        "query_start",
        "query_end",
        "subject_start",
        "subject_end",
        "evalue",
        "bit_score",
    ]
    writer = csv.DictWriter(
        handle,
        fieldnames=fieldnames,
        delimiter="\t",
        lineterminator="\n",
    )
    writer.writeheader()
    for rank, hit in enumerate(hits, start=1):
        writer.writerow(
            {
                "rank": rank,
                "te_id": hit.subject_id,
                "percent_identity": f"{hit.percent_identity:.3f}",
                "query_coverage_pct": f"{hit.query_coverage_pct:.3f}",
                "alignment_length": hit.alignment_length,
                "query_length": hit.query_length,
                "subject_length": hit.subject_length,
                "mismatches": hit.mismatches,
                "gap_opens": hit.gap_opens,
                "query_start": hit.query_start,
                "query_end": hit.query_end,
                "subject_start": hit.subject_start,
                "subject_end": hit.subject_end,
                "evalue": hit.evalue,
                "bit_score": f"{hit.bit_score:.1f}",
            }
        )


def load_input_sequence(args: argparse.Namespace) -> str:
    if args.sequence is not None:
        return sequence_from_text(args.sequence)
    if args.sequence_file is not None:
        return load_sequence_from_file(args.sequence_file.expanduser().resolve())
    return sequence_from_text(sys.stdin.read())


def main(argv: Sequence[str]) -> int:
    args = parse_args(argv)
    te_library = require_file(args.te_library, "TE library FASTA")
    blastn_exe = require_executable(args.blastn, "blastn")
    makeblastdb_exe = require_executable(args.makeblastdb, "makeblastdb")
    seq = load_input_sequence(args)
    db_prefix = build_blast_db(
        te_library=te_library,
        db_dir=args.db_dir,
        makeblastdb_exe=makeblastdb_exe,
    )
    hits = run_blast(
        seq=seq,
        db_prefix=db_prefix,
        blastn_exe=blastn_exe,
        task=args.task,
        evalue=args.evalue,
        top=args.top,
    )
    if args.output is None:
        write_hits_tsv(sys.stdout, hits)
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        with args.output.open("w", encoding="utf-8", newline="") as handle:
            write_hits_tsv(handle, hits)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
