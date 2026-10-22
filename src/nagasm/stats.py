"""Assembly statistics.

N50 and its relatives, worked out here rather than pulled in from another tool.
They are a dozen lines each, and having them in the repo means the numbers in the
report are ones you can read the code for.

A warning about N50 that is worth repeating, because it bites on this organism in
particular. N50 rewards long contigs and knows nothing about whether they are
right. An assembler that glues two chromosomes together scores better than one
that leaves them apart, and an assembly that quietly collapsed the thousands of
copies of the rDNA circle into one long tandem array would look excellent by this
measure and be badly wrong. Read these numbers next to the compartment table, not
instead of it.
"""

from __future__ import annotations

from dataclasses import dataclass, asdict
from typing import Iterable, Sequence

from .seqio import Record


@dataclass
class AssemblyStats:
    n_contigs: int
    total_bp: int
    longest_bp: int
    shortest_bp: int
    mean_bp: float
    median_bp: float
    n50: int
    l50: int
    n90: int
    l90: int
    gc_percent: float
    n_percent: float
    n_circular: int = 0

    def as_dict(self) -> dict:
        return asdict(self)

    def summary_lines(self) -> list[str]:
        """The stats as plain sentences, for a report or a log."""
        return [
            f"Contigs: {self.n_contigs:,}",
            f"Total length: {self.total_bp:,} bp",
            f"Longest contig: {self.longest_bp:,} bp",
            f"N50: {self.n50:,} bp (L50 {self.l50:,})",
            f"N90: {self.n90:,} bp (L90 {self.l90:,})",
            f"GC: {self.gc_percent:.1f}%",
            f"Unknown bases: {self.n_percent:.3f}%",
        ]


def nx(lengths: Sequence[int], fraction: float) -> tuple[int, int]:
    """The Nx length and the Lx count for a set of contig lengths.

    Sort the contigs longest first and add their lengths up. Nx is the length of
    the contig that takes the running total past the given fraction of the whole
    assembly, and Lx is how many contigs that took. So N50 of 700 kb means half
    the assembly sits in contigs of at least 700 kb.
    """
    if not lengths or fraction <= 0 or fraction > 1:
        return 0, 0
    ordered = sorted(lengths, reverse=True)
    target = sum(ordered) * fraction
    running = 0
    for i, length in enumerate(ordered, start=1):
        running += length
        if running >= target:
            return length, i
    return ordered[-1], len(ordered)


def gc_and_n(sequences: Iterable[str]) -> tuple[float, float]:
    """GC share of the readable bases, and the share of bases that are unknown.

    These are two different denominators on purpose. GC is of the bases we could
    read, because an N is not evidence either way. The N figure is of everything,
    because that is the question being asked.
    """
    gc = 0
    acgt = 0
    unknown = 0
    total = 0
    for seq in sequences:
        for base in seq:
            total += 1
            upper = base.upper()
            if upper in "GC":
                gc += 1
                acgt += 1
            elif upper in "AT":
                acgt += 1
            else:
                unknown += 1
    gc_pct = 100.0 * gc / acgt if acgt else 0.0
    n_pct = 100.0 * unknown / total if total else 0.0
    return gc_pct, n_pct


def summarise(records: Sequence[Record], n_circular: int = 0) -> AssemblyStats:
    """Describe an assembly."""
    lengths = [len(r.seq) for r in records]
    if not lengths:
        return AssemblyStats(0, 0, 0, 0, 0.0, 0.0, 0, 0, 0, 0, 0.0, 0.0, 0)

    ordered = sorted(lengths)
    middle = len(ordered) // 2
    median = (
        float(ordered[middle])
        if len(ordered) % 2
        else (ordered[middle - 1] + ordered[middle]) / 2.0
    )
    n50, l50 = nx(lengths, 0.5)
    n90, l90 = nx(lengths, 0.9)
    gc_pct, n_pct = gc_and_n(r.seq for r in records)

    return AssemblyStats(
        n_contigs=len(lengths),
        total_bp=sum(lengths),
        longest_bp=max(lengths),
        shortest_bp=min(lengths),
        mean_bp=sum(lengths) / len(lengths),
        median_bp=median,
        n50=n50,
        l50=l50,
        n90=n90,
        l90=l90,
        gc_percent=gc_pct,
        n_percent=n_pct,
        n_circular=n_circular,
    )


def read_length_summary(reads: Sequence[str]) -> dict:
    """Describe a read set, including read N50.

    Read N50 is the number worth watching for a long read library. The mean is
    dragged down by a crowd of short reads that contribute almost nothing, while
    N50 says how long the reads carrying most of the sequence actually are, and
    that is what decides which repeats can be crossed.
    """
    lengths = [len(r) for r in reads]
    if not lengths:
        return {
            "n_reads": 0,
            "total_bp": 0,
            "mean_bp": 0.0,
            "median_bp": 0.0,
            "n50_bp": 0,
            "longest_bp": 0,
        }
    ordered = sorted(lengths)
    middle = len(ordered) // 2
    median = (
        float(ordered[middle])
        if len(ordered) % 2
        else (ordered[middle - 1] + ordered[middle]) / 2.0
    )
    n50, _ = nx(lengths, 0.5)
    return {
        "n_reads": len(lengths),
        "total_bp": sum(lengths),
        "mean_bp": sum(lengths) / len(lengths),
        "median_bp": median,
        "n50_bp": n50,
        "longest_bp": max(lengths),
    }


def coverage_estimate(total_read_bp: int, genome_bp: int) -> float:
    """How many times over the reads cover a genome of a given size."""
    return total_read_bp / genome_bp if genome_bp else 0.0
