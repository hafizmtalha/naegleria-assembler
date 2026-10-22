"""An overlap assembler of our own, for the nuclear compartment.

This is stage two, and it is honest about what it is. Flye and SPAdes are the
result of many years of work on the hard parts of assembly, and this is not going
to beat them. It is here so that someone can read an assembler end to end in a
few hundred lines and see what the steps actually are, and so that the package
does not depend on an external tool to produce an answer at all.

The approach is overlap, layout, consensus, which is the oldest way of doing this
and the easiest to follow:

1. Find which reads overlap, using the anchors from the Rust core.
2. Throw away overlaps that look like they come from repeats rather than from two
   reads of the same place.
3. Walk chains of reads, each one extending the last, until nothing extends.
4. Write out each chain as a contig.

Compared with the circular assembler next door, the difference is that a
chromosome has ends. The walk here runs out of reads and stops, which is the
normal outcome, and then starts again somewhere else. Every read is used at most
once, so the contigs do not overlap each other.

## Where this falls short, said plainly

Repeats are the hard part of assembly and this handles them bluntly. An overlap
whose anchors are spread thinly, or a read that overlaps a crowd of others at the
same place, is dropped rather than resolved. That breaks the assembly at every
repeat longer than a read instead of working out which copy is which. The result
is more contigs than a good assembler would give, but the contigs it does give
are not glued together wrongly, which is the failure that matters more.
"""

from __future__ import annotations

from dataclasses import dataclass

from . import nagcore
from .circular import _overlap_graph, _place_neighbour
from .seqio import Record, reverse_complement


@dataclass
class Contig:
    sequence: str
    read_indices: list[int]

    def __len__(self) -> int:
        return len(self.sequence)


def _drop_crowded_reads(
    graph: dict[int, list[tuple[int, bool, int, int]]], max_neighbours: int
) -> set[int]:
    """Reads that overlap far more partners than coverage can explain.

    A read sitting inside a repeat matches every copy of that repeat, so it has a
    crowd of partners that do not belong together. Chaining through it welds
    unrelated parts of the genome into one contig, which is the worst thing an
    assembler can do, because nothing downstream will tell you it happened.
    Leaving these reads out breaks the assembly at repeats instead, which is
    visible and safe.
    """
    return {read for read, neighbours in graph.items() if len(neighbours) > max_neighbours}


def _walk_chain(
    seed: int,
    reads: list[str],
    graph: dict[int, list[tuple[int, bool, int, int]]],
    rightwards: bool,
    used: set[int],
    blocked: set[int],
    k: int,
) -> list[tuple[int, int, bool]]:
    """Extend from a seed in one direction. Returns (read, offset, forward)."""
    out: list[tuple[int, int, bool]] = []
    current, offset, forward = seed, 0, True

    while True:
        best: tuple[int, int, bool] | None = None
        best_reach: int | None = None
        for other, same_strand, shift, _anchors in graph.get(current, ()):
            if other in used or other in blocked:
                continue
            new_offset, new_forward = _place_neighbour(
                offset, forward, len(reads[current]), len(reads[other]), same_strand, shift, k
            )
            advance = new_offset - offset
            if not rightwards:
                advance = -advance
            if advance <= 0 or advance >= len(reads[current]):
                continue
            reach = new_offset if rightwards else -new_offset
            if best_reach is None or reach > best_reach:
                best, best_reach = (other, new_offset, new_forward), reach

        if best is None:
            return out
        other, new_offset, new_forward = best
        used.add(other)
        out.append(best)
        current, offset, forward = other, new_offset, new_forward


def _stitch(reads: list[str], placements: list[tuple[int, int, bool]]) -> str:
    """Lay the chain's reads end to end, each adding what the last did not cover."""
    if not placements:
        return ""
    ordered = sorted(placements, key=lambda p: p[1])
    base = ordered[0][1]
    pieces: list[str] = []
    covered = 0
    for read_id, offset, forward in ordered:
        seq = reads[read_id]
        if not forward:
            seq = reverse_complement(seq)
        start = offset - base
        end = start + len(seq)
        if end <= covered:
            continue
        pieces.append(seq[max(0, covered - start) :])
        covered = end
    return "".join(pieces)


def assemble_olc(
    reads: list[str],
    k: int = 15,
    w: int = 10,
    min_entropy: float = 1.0,
    min_anchors: int = 5,
    min_shared: int = 4,
    max_neighbours: int = 200,
    min_contig_bp: int = 500,
) -> list[Record]:
    """Assemble reads into contigs with overlap, layout and consensus.

    The draft contigs are raw read sequence joined at the seams, so on noisy long
    reads they carry the error rate of the reads. Polish them afterwards, the
    same as for the circular assembler.
    """
    reads = [r for r in reads if len(r) >= k * 2]
    if not reads:
        return []

    graph = _overlap_graph(reads, k, w, min_entropy, min_anchors, min_shared)
    if not graph:
        # Nothing overlaps anything, so every read is its own contig.
        return [
            Record(f"contig_{i + 1}", r)
            for i, r in enumerate(sorted(reads, key=len, reverse=True))
            if len(r) >= min_contig_bp
        ]

    blocked = _drop_crowded_reads(graph, max_neighbours)
    used: set[int] = set()
    contigs: list[Contig] = []

    # Start each chain from the read with the most overlap support that has not
    # been used yet, so the strongest parts of the graph are laid out first and
    # the leftovers fill in around them.
    order = sorted(
        (r for r in graph if r not in blocked),
        key=lambda r: sum(n[3] for n in graph[r]),
        reverse=True,
    )

    for seed in order:
        if seed in used:
            continue
        used.add(seed)
        chain = [(seed, 0, True)]
        chain += _walk_chain(seed, reads, graph, True, used, blocked, k)
        chain += _walk_chain(seed, reads, graph, False, used, blocked, k)
        sequence = _stitch(reads, chain)
        if sequence:
            contigs.append(Contig(sequence, [c[0] for c in chain]))

        # Every read overlapping this chain is already inside the sequence we
        # just wrote, so retire them all.
        #
        # Without this the assembler produces the genome over and over. At
        # thirty fold coverage only about one read in thirty ends up on the
        # chain, and each of the others is still free to seed a chain of its own
        # that runs along the same stretch of genome. An early version of this
        # function turned a 120 kb genome into 2.8 Mb spread over a hundred
        # contigs, every one of them a correct and nearly complete copy of the
        # same thing.
        for read_id, _, _ in chain:
            for other, *_ in graph.get(read_id, ()):
                used.add(other)

    # Reads that were blocked or never joined anything still hold sequence.
    for i, read in enumerate(reads):
        if i not in used and len(read) >= min_contig_bp:
            contigs.append(Contig(read, [i]))

    contigs = [c for c in contigs if len(c) >= min_contig_bp]
    contigs.sort(key=len, reverse=True)
    return [Record(f"contig_{i + 1}", c.sequence) for i, c in enumerate(contigs)]


def overlap_summary(reads: list[str], k: int = 15, w: int = 10) -> dict:
    """A quick look at the overlap graph, without assembling anything.

    Useful for working out whether coverage is good enough before spending time
    on an assembly. A graph where most reads have no partners means the reads are
    too few, too short or too noisy for the k being used.
    """
    index = nagcore.OverlapIndex(k=k, w=w, min_entropy=1.0)
    index.add_reads(reads)
    hits = index.place_overlaps(min_shared=3, max_occurrence=10_000)
    partners: dict[int, int] = {}
    for a, b, *_ in hits:
        partners[a] = partners.get(a, 0) + 1
        partners[b] = partners.get(b, 0) + 1
    counts = [partners.get(i, 0) for i in range(len(reads))]
    isolated = sum(1 for c in counts if c == 0)
    return {
        "n_reads": len(reads),
        "n_overlaps": len(hits),
        "reads_with_no_overlap": isolated,
        "share_isolated": isolated / len(reads) if reads else 0.0,
        "mean_partners": sum(counts) / len(counts) if counts else 0.0,
        "max_partners": max(counts) if counts else 0,
    }
