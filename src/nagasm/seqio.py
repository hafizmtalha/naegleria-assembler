"""Reading and writing FASTA and FASTQ, including gzipped files.

Nothing clever here. A small reader keeps the rest of the package free of a
dependency that would otherwise show up in every module, and it means the tests
can make a file and read it back without installing anything.
"""

from __future__ import annotations

import gzip
import io
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Iterator


@dataclass(slots=True)
class Record:
    """One sequence, with its name and, for FASTQ, its quality string."""

    name: str
    seq: str
    qual: str | None = None

    def __len__(self) -> int:
        return len(self.seq)


def _open_text(path: str | Path) -> io.TextIOBase:
    """Open a file, transparently handling gzip."""
    path = Path(path)
    with open(path, "rb") as probe:
        magic = probe.read(2)
    if magic == b"\x1f\x8b":
        return gzip.open(path, "rt")
    return open(path, "rt")


def read_fasta(path: str | Path) -> Iterator[Record]:
    """Yield every record in a FASTA file.

    Sequence lines are joined, so wrapped files read the same as unwrapped ones.
    """
    name: str | None = None
    chunks: list[str] = []
    with _open_text(path) as handle:
        for line in handle:
            line = line.rstrip("\n\r")
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    yield Record(name, "".join(chunks))
                name = line[1:].split()[0] if len(line) > 1 else ""
                chunks = []
            else:
                chunks.append(line)
    if name is not None:
        yield Record(name, "".join(chunks))


def read_fastq(path: str | Path) -> Iterator[Record]:
    """Yield every record in a FASTQ file.

    This expects the usual four lines per record. Sequences wrapped over several
    lines are rare enough in real FASTQ that supporting them would cost more
    clarity than it is worth, and a malformed file raises rather than guessing.
    """
    with _open_text(path) as handle:
        while True:
            header = handle.readline()
            if not header:
                return
            header = header.rstrip("\n\r")
            if not header.startswith("@"):
                raise ValueError(f"expected a line starting with @, found: {header[:60]!r}")
            seq = handle.readline().rstrip("\n\r")
            plus = handle.readline()
            if not plus.startswith("+"):
                raise ValueError(f"expected a line starting with +, found: {plus[:60]!r}")
            qual = handle.readline().rstrip("\n\r")
            if len(qual) != len(seq):
                raise ValueError(
                    f"{header[1:].split()[0]}: {len(seq)} bases but {len(qual)} quality values"
                )
            yield Record(header[1:].split()[0], seq, qual)


def read_sequences(path: str | Path) -> Iterator[Record]:
    """Read a file whether it is FASTA or FASTQ, picking by its first character."""
    with _open_text(path) as handle:
        first = handle.read(1)
    if first == ">":
        yield from read_fasta(path)
    elif first == "@":
        yield from read_fastq(path)
    elif first == "":
        return
    else:
        raise ValueError(f"{path} does not start with > or @, so it is neither FASTA nor FASTQ")


def write_fasta(path: str | Path, records: Iterable[Record], wrap: int = 60) -> int:
    """Write records to FASTA and give back how many were written.

    Lines are wrapped at `wrap` characters, or left whole if `wrap` is 0.
    """
    n = 0
    with open(path, "wt") as handle:
        for rec in records:
            handle.write(f">{rec.name}\n")
            seq = rec.seq
            if wrap and wrap > 0:
                for i in range(0, len(seq), wrap):
                    handle.write(seq[i : i + wrap] + "\n")
            else:
                handle.write(seq + "\n")
            n += 1
    return n


def write_fastq(path: str | Path, records: Iterable[Record]) -> int:
    """Write records to FASTQ. A record with no quality gets a flat placeholder."""
    n = 0
    with open(path, "wt") as handle:
        for rec in records:
            qual = rec.qual if rec.qual is not None else "I" * len(rec.seq)
            handle.write(f"@{rec.name}\n{rec.seq}\n+\n{qual}\n")
            n += 1
    return n


_COMPLEMENT = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")


def reverse_complement(seq: str) -> str:
    """The other strand, read the usual way round."""
    return seq.translate(_COMPLEMENT)[::-1]
