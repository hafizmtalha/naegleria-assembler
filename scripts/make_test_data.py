#!/usr/bin/env python3
"""Make a small read set you can run every example on.

Real reads are too big to keep in a git repository, and the public ones need
downloading before anything works. So this builds a fake genome instead, writes
reads from it, and the whole thing takes a few seconds and about 20 MB.

The point is not that the sequence is realistic. It is that the **shape** is. The
pipeline exists to tell three things apart, so the fake genome holds one of each:

1. **A chromosome**, long and at ordinary depth. This is what Flye should
   assemble and what everything else gets measured against.

2. **A circular element at high copy.** Far deeper than the chromosome and much
   shorter, like a plasmid or the rDNA circle in Naegleria. The depth sorting
   should pull it out into its own compartment, the circle finder should close
   it, and the support check should accept it. This is the thing the pipeline is
   for.

3. **A repeat family inside the chromosome.** The same block, several times
   over, scattered along it. This one is here to make the chromosome realistic
   rather than to test the support check, and the next section says why.

If the run finds the circle, closes it at about the right length, and the support
check accepts it, the pipeline is working.

## What this data does not test

It does not reproduce the failure the support check exists to catch, and it is
worth being straight about that rather than implying otherwise.

On real E. coli data the pipeline reports several elements that are not real, and
they score 0.10 to 0.14 on containment, well under the 0.20 threshold, so they
get rejected. Those are **mosaics**: the overlap walk stitched together fragments
from several different copies of a repeat, and the result is a sequence that no
single read covers more than about a tenth of. That is what the low score means
and that is what makes them catchable.

A simulated repeat is one identical block. Any read that touches it covers the
whole block, so the share works out as repeat length over read length. At the
1,200 bp default in 9 kb reads that is about 0.13, which looks right, but those
reads never get that far: a read is assigned to a compartment by the median depth
of its k-mers, and 7,800 bp of unique chromosome outvotes 1,200 bp of repeat. So
the reads stay in the nuclear compartment, which is the correct thing to do with
them, and no fake element is ever built.

Make the repeat long enough to dominate a read and it does form a compartment,
but then the share genuinely rises above the threshold and it passes, because at
that point a third of every read really is inside it. Try
`--repeat 3000 --repeat-copies 8` to see that happen: the check reports 0.32 and
accepts it, which is honest arithmetic rather than a bug.

Building a mosaic on purpose would mean several similar but distinct repeat
families and would depend on how the assembler happens to walk them, which is not
something a generator can promise. So the rejection path is covered by the real
worked example in the README, not here.

Reads carry errors, because clean reads hide most of the problems worth finding.
The default of 6 percent is roughly nanopore, with indels running about twice as
often as substitutions, which is the part that actually breaks assemblers.

Usage:

    python scripts/make_test_data.py                     # into data/simulated
    python scripts/make_test_data.py --out somewhere/else
    python scripts/make_test_data.py --error 0.12        # harder
    python scripts/make_test_data.py --small             # quicker, no Flye needed

Everything is seeded, so the same arguments give the same reads every time. Pass
a different --seed if you want a second, different set.
"""

from __future__ import annotations

import argparse
import gzip
import random
from pathlib import Path

BASES = "ACGT"


def random_sequence(n: int, rng: random.Random, gc: float = 0.36) -> str:
    """A sequence of n bases at roughly the GC fraction asked for.

    The default of 0.36 is close to Naegleria, which is AT rich. It matters a
    little here because the depth sorting looks at GC, and a genome sitting at
    exactly 50 percent would not show the spread a real one does.
    """
    out = []
    for _ in range(n):
        out.append(rng.choice("GC") if rng.random() < gc else rng.choice("AT"))
    return "".join(out)


def with_repeats(
    backbone: str, repeat: str, copies: int, rng: random.Random
) -> tuple[str, list[int]]:
    """Drop several copies of one block into a chromosome, spaced out.

    Returns the new sequence and where the copies landed, so the caller can say
    in the notes where they are. Spacing them out matters. Copies sitting next to
    each other make a tandem array, which is a different problem, and what we
    want here is the scattered kind that a naive assembler stitches together
    into a fake circle.
    """
    n = len(backbone)
    step = n // (copies + 1)
    pieces = []
    at = []
    cursor = 0
    for i in range(copies):
        # Put each copy near its own slot, jittered so they are not on a grid.
        target = step * (i + 1) + rng.randint(-step // 6, step // 6)
        target = max(cursor, min(target, n - 1))
        pieces.append(backbone[cursor:target])
        at.append(sum(len(p) for p in pieces))
        pieces.append(repeat)
        cursor = target
    pieces.append(backbone[cursor:])
    return "".join(pieces), at


def mutate(sequence: str, error: float, rng: random.Random) -> str:
    """Put sequencing errors into one read.

    Indels get two thirds of the budget and substitutions one third, because that
    is the balance on nanopore data and it is the indels that cause the trouble.
    A substitution is a wrong letter in the right place, which assemblers shrug
    off. An indel shifts everything after it.
    """
    out = []
    for base in sequence:
        roll = rng.random()
        if roll >= error:
            out.append(base)
            continue
        kind = rng.random()
        if kind < 1 / 3:
            out.append(rng.choice(BASES.replace(base, "") if base in BASES else BASES))
        elif kind < 2 / 3:
            pass  # a deletion
        else:
            out.append(base)
            out.append(rng.choice(BASES))  # an insertion
    return "".join(out)


def sample_reads(
    sequence: str,
    depth: float,
    mean_length: int,
    error: float,
    rng: random.Random,
    circular: bool = False,
) -> list[str]:
    """Reads from one molecule, at roughly the depth asked for.

    For a circle the sequence is doubled before sampling, so reads can run across
    the join. Without that every read stops at the same place, the join never gets
    covered, and the circle finder has nothing to work with. That is an artefact
    of the simulation rather than anything real, so it has to be avoided here.
    """
    source = sequence + sequence if circular else sequence
    span = len(sequence)
    wanted = int(span * depth / mean_length)
    reads = []
    for _ in range(max(wanted, 1)):
        # Lengths vary a lot on long read data, so this spreads them widely
        # rather than clustering at the mean.
        length = max(500, int(rng.gauss(mean_length, mean_length * 0.45)))
        length = min(length, span if not circular else span * 2 - 1)
        start = rng.randrange(0, span if circular else max(span - length, 1))
        piece = source[start : start + length]
        if len(piece) < 200:
            continue
        if rng.random() < 0.5:
            piece = reverse_complement(piece)
        reads.append(mutate(piece, error, rng))
    return reads


def reverse_complement(sequence: str) -> str:
    table = str.maketrans("ACGTN", "TGCAN")
    return sequence.translate(table)[::-1]


def write_fastq(path: Path, reads: list[str], prefix: str, quality: str = "I") -> None:
    """Write reads as gzipped FASTQ with flat quality.

    The quality string is the same letter all the way along. Nothing in this
    pipeline reads quality scores, so inventing a believable pattern would be
    effort spent on something no code looks at.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wt") as handle:
        for i, read in enumerate(reads, 1):
            handle.write(f"@{prefix}_{i} length={len(read)}\n{read}\n+\n{quality * len(read)}\n")


def write_fasta(path: Path, records: list[tuple[str, str]], wrap: int = 70) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as handle:
        for name, sequence in records:
            handle.write(f">{name}\n")
            if wrap:
                for i in range(0, len(sequence), wrap):
                    handle.write(sequence[i : i + wrap] + "\n")
            else:
                handle.write(sequence + "\n")


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Make a small read set with a high copy circle in it, so the "
        "examples in the README have something to run on.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "-o", "--out", default="data/simulated", help="where to write (default: data/simulated)"
    )
    parser.add_argument(
        "--chromosome", type=int, default=400_000, help="chromosome length in bp (default: 400000)"
    )
    parser.add_argument(
        "--element", type=int, default=6_000, help="circle length in bp (default: 6000)"
    )
    parser.add_argument(
        "--repeat",
        type=int,
        default=1_200,
        help="repeat block length in bp (default: 1200, about the size of an IS "
        "element). Raising this past roughly a fifth of --read-length changes "
        "what the data demonstrates, so read the module docstring first",
    )
    parser.add_argument(
        "--repeat-copies", type=int, default=8, help="copies of the repeat (default: 8)"
    )
    parser.add_argument(
        "--depth", type=float, default=35.0, help="chromosome depth (default: 35)"
    )
    parser.add_argument(
        "--element-depth",
        type=float,
        default=900.0,
        help="circle depth. Must be well above the chromosome or the depth sorting "
        "has nothing to separate (default: 900)",
    )
    parser.add_argument(
        "--read-length", type=int, default=9_000, help="mean read length (default: 9000)"
    )
    parser.add_argument(
        "--error", type=float, default=0.06, help="per base error rate (default: 0.06)"
    )
    parser.add_argument("--seed", type=int, default=1, help="for repeatable output (default: 1)")
    parser.add_argument(
        "--small",
        action="store_true",
        help="a quarter the size, for a quick check. Too small for Flye, so use "
        "--skip-nuclear with it",
    )
    args = parser.parse_args()

    if args.small:
        args.chromosome //= 4
        args.depth = 25.0
        args.element_depth = 600.0

    if args.element_depth <= args.depth * 4:
        parser.error(
            f"the circle depth ({args.element_depth:.0f}) has to be well above the "
            f"chromosome depth ({args.depth:.0f}), or there is no hump for the depth "
            "sorting to find and this data will not demonstrate anything"
        )

    rng = random.Random(args.seed)
    out = Path(args.out)

    # The repeat block, and a chromosome with several copies of it inside.
    repeat = random_sequence(args.repeat, rng)
    backbone = random_sequence(args.chromosome, rng)
    chromosome, repeat_at = with_repeats(backbone, repeat, args.repeat_copies, rng)

    # The circle. Deliberately unlike the chromosome, so a match to it means
    # something.
    element = random_sequence(args.element, rng, gc=0.44)

    print("Building:")
    print(f"  chromosome   {len(chromosome):>9,} bp at {args.depth:.0f}x")
    print(f"               with {args.repeat_copies} copies of a {args.repeat:,} bp repeat")
    print(f"  circle       {len(element):>9,} bp at {args.element_depth:.0f}x")
    print(f"  read length  {args.read_length:>9,} bp mean, {args.error:.0%} error")

    reads = sample_reads(chromosome, args.depth, args.read_length, args.error, rng)
    reads += sample_reads(
        element, args.element_depth, args.read_length, args.error, rng, circular=True
    )
    rng.shuffle(reads)

    reads_path = out / "reads.fastq.gz"
    write_fastq(reads_path, reads, "sim")

    # The answers, so anyone can check what the run reported against what went in.
    write_fasta(
        out / "truth_chromosome.fasta",
        [("chromosome", chromosome)],
    )
    write_fasta(out / "truth_element.fasta", [("element_circular", element)])
    write_fasta(out / "truth_repeat.fasta", [("repeat_block", repeat)])

    total = sum(len(r) for r in reads)
    (out / "README.md").write_text(
        f"""# Simulated test data

Made by `scripts/make_test_data.py` with seed {args.seed}. Delete it and run the
script again to get it back, byte for byte.

| | |
|---|---|
| Reads | {len(reads):,} |
| Bases | {total:,} |
| Chromosome | {len(chromosome):,} bp at {args.depth:.0f}x |
| Repeat inside it | {args.repeat:,} bp, {args.repeat_copies} copies |
| Circle | {len(element):,} bp at {args.element_depth:.0f}x |
| Error rate | {args.error:.0%} |

## What a good run looks like

The circle is the thing to watch. A run that works finds a high copy compartment,
closes a circle of about {len(element):,} bp from it, and the support check
accepts it with a containment score near 0.99. The chromosome should come back at
close to {len(chromosome):,} bp.

The repeat inside the chromosome is there to make the assembly realistic, not to
test the support check. Reads that touch it stay in the nuclear compartment,
because a read is placed by the median depth of its k-mers and the unique
sequence either side of the repeat outvotes it. That is the right thing to do
with those reads, so no fake element gets built and there is nothing for the
check to reject.

If you want to watch the support check work, use the real worked example in the
main README instead. On E. coli it rejects four elements at 0.10 to 0.14
containment. Those are mosaics of several repeat copies, which a generator cannot
promise to recreate, and `scripts/make_test_data.py` explains why in more detail.

## Files

| File | What it is |
|---|---|
| `reads.fastq.gz` | the reads, which is all the pipeline needs |
| `truth_chromosome.fasta` | the chromosome the reads came from |
| `truth_element.fasta` | the circle, to compare what you assembled against |
| `truth_repeat.fasta` | the repeat block, to see where it turns up |

The truth files are for checking afterwards. Do not pass them to `--reference`
while you are testing whether the pipeline can find things on its own, since that
is giving it the answer.
"""
    )

    print()
    print(f"Wrote {len(reads):,} reads, {total:,} bases to {reads_path}")
    print(f"  repeat copies sit at: {', '.join(f'{p:,}' for p in repeat_at)}")
    print()
    print("Try:")
    print(f"  nagasm triage {reads_path} --scale 20")
    print(f"  nagasm assemble -l {reads_path} -o results/simulated -t 8 --scale 20")
    if args.small:
        print("  (made with --small, so add --skip-nuclear, it is too small for Flye)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
