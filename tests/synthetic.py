"""Making read sets with a known answer, so the tests can check the answer.

A real Naegleria read set has a shape that is hard to reason about from the
outside: a few tens of megabases of chromosome at ordinary depth, with two small
circles buried in it at hundreds or thousands of times that depth. Everything in
this package is built around finding that shape, so the tests need to be able to
produce it on demand at a size that runs in a second.

Nothing here is a model of Naegleria biology. It is a way of making a pile of
reads whose compartments we already know, so a test can ask whether the code
found them.
"""

from __future__ import annotations

import random
from dataclasses import dataclass

_BASES = "ACGT"
_COMPLEMENT = str.maketrans("ACGT", "TGCA")


def random_sequence(length: int, rng: random.Random, gc: float = 0.5) -> str:
    """Random DNA of a given length at roughly a given GC share.

    GC matters here because Naegleria is AT rich, and several parts of the
    package behave differently on AT rich sequence. A test that only ever used an
    even base mix would miss that.
    """
    weights = [(1 - gc) / 2, gc / 2, gc / 2, (1 - gc) / 2]
    return "".join(rng.choices(_BASES, weights=weights, k=length))


def reverse_complement(seq: str) -> str:
    return seq.translate(_COMPLEMENT)[::-1]


def add_errors(seq: str, rate: float, rng: random.Random) -> str:
    """Scatter substitutions, insertions and deletions through a sequence.

    Long read platforms make all three, in roughly similar amounts, and the
    homopolymer handling in this package exists because of the indels. Errors
    here are split evenly between the three kinds.
    """
    if rate <= 0:
        return seq
    out: list[str] = []
    for base in seq:
        roll = rng.random()
        if roll >= rate:
            out.append(base)
        elif roll < rate / 3:
            out.append(rng.choice(_BASES))  # substitution
        elif roll < 2 * rate / 3:
            out.append(base)
            out.append(rng.choice(_BASES))  # insertion
        # the remaining third is a deletion, so nothing is appended
    return "".join(out)


def sample_reads(
    template: str,
    depth: float,
    mean_length: int,
    rng: random.Random,
    circular: bool = False,
    error_rate: float = 0.0,
    min_length: int = 200,
) -> list[str]:
    """Cut a template into reads at roughly the depth asked for.

    When `circular` is true, reads are allowed to run off the end and wrap round
    to the start, which is what actually happens on a circular molecule and what
    makes closing a circle possible at all.
    """
    if not template:
        return []
    n_reads = max(1, int(round(depth * len(template) / mean_length)))
    reads: list[str] = []
    for _ in range(n_reads):
        length = max(min_length, int(rng.gauss(mean_length, mean_length * 0.3)))
        if circular:
            start = rng.randrange(len(template))
            doubled = template + template
            length = min(length, len(template))
            piece = doubled[start : start + length]
        else:
            if length >= len(template):
                piece = template
            else:
                start = rng.randrange(0, len(template) - length + 1)
                piece = template[start : start + length]
        if rng.random() < 0.5:
            piece = reverse_complement(piece)
        if error_rate:
            piece = add_errors(piece, error_rate, rng)
        if len(piece) >= min_length:
            reads.append(piece)
    return reads


@dataclass
class SyntheticGenome:
    """A read set plus the truth about where every read came from."""

    reads: list[str]
    labels: list[str]
    templates: dict[str, str]
    depths: dict[str, float]

    def reads_for(self, label: str) -> list[str]:
        return [r for r, lab in zip(self.reads, self.labels) if lab == label]


def naegleria_like(
    seed: int = 0,
    nuclear_bp: int = 400_000,
    nuclear_depth: float = 25,
    mito_bp: int = 12_000,
    mito_depth: float = 300,
    rdna_bp: int = 4_000,
    rdna_depth: float = 2_000,
    mean_read_length: int = 2_000,
    error_rate: float = 0.0,
    gc: float = 0.36,
) -> SyntheticGenome:
    """A small read set with the shape of a Naegleria library.

    The sizes are scaled right down so tests run quickly. What is kept is the
    thing that matters: a large compartment at low depth and two small ones at far
    higher depth, with the deepest compartment contributing a lot of reads
    despite being tiny.

    The GC default is near the published figure for Naegleria, so the AT rich
    behaviour gets exercised.
    """
    rng = random.Random(seed)
    templates = {
        "nuclear": random_sequence(nuclear_bp, rng, gc=gc),
        "mito": random_sequence(mito_bp, rng, gc=gc - 0.03),
        "rdna": random_sequence(rdna_bp, rng, gc=gc + 0.08),
    }
    depths = {"nuclear": nuclear_depth, "mito": mito_depth, "rdna": rdna_depth}
    circular = {"nuclear": False, "mito": True, "rdna": True}

    reads: list[str] = []
    labels: list[str] = []
    for name, template in templates.items():
        got = sample_reads(
            template,
            depths[name],
            mean_read_length,
            rng,
            circular=circular[name],
            error_rate=error_rate,
        )
        reads.extend(got)
        labels.extend([name] * len(got))

    order = list(range(len(reads)))
    rng.shuffle(order)
    return SyntheticGenome(
        reads=[reads[i] for i in order],
        labels=[labels[i] for i in order],
        templates=templates,
        depths=depths,
    )
