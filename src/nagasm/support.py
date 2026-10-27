"""Asking whether an assembled element is really a separate molecule.

The depth sorting finds sequence that is present in more copies than the main
genome. In Naegleria that works, because the deep sequence really is held on
separate small circles. The step that turns a deep compartment into a circle
quietly assumes the same thing everywhere, and that assumption is wrong often
enough to matter.

A bacterium is where it breaks. Klebsiella carries eight copies of its
ribosomal RNA operon and a scatter of insertion sequences, all inside the
chromosome. That sequence is genuinely several times deeper than the rest of
the genome, so the depth sorting puts it in a high copy compartment and is
right to. It is simply not a separate molecule. Handed those reads, the
circular assembler closes loops that do not exist, and we saw it report four of
them on real Klebsiella data, one of which was a 6.5 kb sequence stitched from
reads belonging to dozens of unrelated places in the chromosome.

Nothing in the package caught that, because nothing asked the obvious question
afterwards: do the reads actually support this?

## Two numbers, because one is not enough

**Depth.** Map the reads back and see how deeply the element is covered. An
element pulled out of a compartment sitting at hundreds of times the main
genome cannot honestly come back covered twice. Three of the four Klebsiella
circles were supported at under 3x and were caught here.

Depth on its own is nowhere near sufficient. On real E. coli data every
element came back between 415x and 5,420x, the two genuine plasmids and the
four artefacts alike, so no depth cut separates them.

**Containment.** For each read, what share of it lies inside the element. A
read from a free circle matches along its whole length, because there is
nothing else on that molecule to match. A read crossing the boundary of a
chromosomal repeat matches only in the middle, with unique chromosome either
side that has nowhere to go. Confirmed elements score 0.41 to 0.93. Confirmed
artefacts score 0.00 to 0.14.

Measuring containment needs two corrections or it reports nonsense:

* A circular element has to be doubled first. A read spanning the join aligns
  as two fragments against a circle written out as a line, and would look like
  a read that half matches.
* The share has to be added up across all of a read's alignments, not taken
  from the longest one. A read off a small circle wraps around it several
  times and comes back as several pieces that between them cover nearly all of
  it. Taking the longest piece alone would mark every small element as poorly
  contained no matter how real it is.

Dividing by the full read length is the part that does the work, and the
obvious looking alternative does not. Normalising by whichever is shorter, the
read or the element, lets a 3 kb repeat inside a chromosome score 0.77,
because a 12 kb read matching the repeat end to end counts as fully contained
even though the other 9 kb of it belongs elsewhere. That is exactly the case
this module exists to catch, and a synthetic test built from a chromosome with
eight copies of one repeat is what exposed it.

Neither number is a verdict on its own and the thresholds here are not precise
instruments. They are set to pass everything we have confirmed and reject
everything we have disproved, and an element that fails is reported with the
reason rather than dropped, because a tool that silently deletes findings is
worse than one that hands you a doubtful one with the numbers attached.
"""

from __future__ import annotations

import subprocess
import tempfile
from dataclasses import dataclass, asdict
from pathlib import Path

from .seqio import Record, write_fasta, write_fastq
from .polish import have, _preset_for


@dataclass
class Support:
    """How well the reads back up an element we have assembled."""

    n_reads: int
    depth: float
    contained: float
    measured: bool = True
    note: str = ""

    def as_dict(self) -> dict:
        return asdict(self)


def measure(
    sequence: str,
    reads: list[str],
    is_circular: bool = True,
    read_type: str = "ont",
    threads: int = 4,
    workdir: Path | None = None,
    max_reads: int = 8_000,
) -> Support:
    """Map the reads back to one element and measure depth and containment.

    Returns `measured=False` if minimap2 is not installed, so that a missing
    tool leaves the element unjudged rather than condemned.
    """
    if not sequence or not reads:
        return Support(0, 0.0, 0.0, measured=False, note="nothing to measure")
    if not have("minimap2"):
        return Support(0, 0.0, 0.0, measured=False, note="minimap2 not installed")

    true_length = len(sequence)
    # Doubled, not tiled. Doubling is enough for a read to cross the join.
    #
    # Tiling the circle until it outran the longest read was the obvious way to
    # stop a long read looking half unaligned against a short element, and it
    # worked, but a tiled circle is a perfect tandem repeat and minimap2 takes
    # pathologically long to chain a read against one. It ran for eleven hours
    # on a single Naegleria read set before being stopped.
    #
    # The ceiling it was there to remove is handled in the arithmetic instead,
    # by adding up how much of each read is covered across all of its
    # alignments rather than looking at the longest one on its own. A read off
    # a real circle wraps and comes back as several pieces that between them
    # cover nearly all of it. A read straddling the edge of a chromosomal
    # repeat comes back as one piece covering only the repeat.
    target = sequence + sequence if is_circular else sequence

    # A few thousand reads settle this. Mapping every read of a deep library
    # against every element costs minutes and changes no verdict.
    #
    # Depth has to be scaled back up afterwards. Taking one read in twelve and
    # reporting the depth of the sample would call a 240 fold element a 20 fold
    # one, and the thinness test would start throwing out real elements in
    # exactly the deep libraries this package is built for.
    sampling_factor = 1.0
    if len(reads) > max_reads:
        step = len(reads) / max_reads
        sampling_factor = step
        reads = [reads[int(i * step)] for i in range(max_reads)]

    tmp = Path(workdir) if workdir else Path(tempfile.mkdtemp(prefix="nagasm_support_"))
    tmp.mkdir(parents=True, exist_ok=True)
    target_path = tmp / "element.fasta"
    reads_path = tmp / "reads.fastq"
    write_fasta(target_path, [Record("element", target)], wrap=0)
    write_fastq(reads_path, [Record(f"read{i}", r) for i, r in enumerate(reads)])

    proc = subprocess.run(
        [
            "minimap2",
            "-x",
            _preset_for(read_type),
            "-t",
            str(threads),
            str(target_path),
            str(reads_path),
        ],
        stdout=subprocess.PIPE,
        stderr=subprocess.DEVNULL,
    )
    if proc.returncode != 0:
        return Support(0, 0.0, 0.0, measured=False, note="minimap2 did not run")

    aligned_bp = 0
    # Every stretch of each read that found a home, gathered per read so that
    # the pieces of a wrapping read can be added up.
    pieces: dict[str, list[tuple[int, int]]] = {}
    read_lengths: dict[str, int] = {}
    for line in proc.stdout.decode().splitlines():
        fields = line.split("\t")
        if len(fields) < 11:
            continue
        name = fields[0]
        start, end = int(fields[2]), int(fields[3])
        aligned_bp += end - start
        pieces.setdefault(name, []).append((start, end))
        read_lengths[name] = int(fields[1])

    if not pieces:
        return Support(0, 0.0, 0.0, note="no reads mapped back to this element")

    shares: list[float] = []
    for name, spans in pieces.items():
        read_length = read_lengths[name]
        if read_length <= 0:
            continue
        # Merge the overlapping pieces so that sequence covered twice is not
        # counted twice, which would let a read look more than fully contained.
        covered = 0
        current_start, current_end = -1, -1
        for start, end in sorted(spans):
            if start > current_end:
                covered += max(0, current_end - current_start)
                current_start, current_end = start, end
            else:
                current_end = max(current_end, end)
        covered += max(0, current_end - current_start)
        shares.append(min(1.0, covered / read_length))

    depth = aligned_bp * sampling_factor / true_length
    contained = sum(shares) / len(shares)
    return Support(int(round(len(pieces) * sampling_factor)), depth, contained)


def verdict(
    support: Support,
    min_depth: float = 10.0,
    min_contained: float = 0.20,
) -> tuple[bool, str]:
    """Whether this looks like a real separate molecule, and why.

    `min_depth` is deliberately low. The job here is to throw out elements with
    no support at all, not to insist on deep coverage, because a real element
    in a shallow library is still a real element.

    `min_contained` is the number that catches sequence stitched together from
    scattered copies of a repeat, which is the failure that depth alone misses.
    """
    if not support.measured:
        return True, support.note or "support was not measured"

    if support.n_reads == 0:
        return False, "no reads map back to it, so there is nothing supporting it"

    if support.depth < min_depth:
        return False, (
            f"only {support.depth:.1f} fold covered by the reads it came from, which "
            f"is too thin to be a real element"
        )

    if support.contained < min_contained:
        return False, (
            f"the reads supporting it only lie {support.contained:.0%} inside it, so "
            f"most of each read belongs somewhere else. This is what sequence "
            f"stitched together from scattered copies of a repeat looks like, rather "
            f"than a separate molecule"
        )

    return True, (
        f"{support.depth:.0f} fold covered, and its reads lie {support.contained:.0%} "
        f"inside it"
    )
