"""Tests for the overlap assembler in this package.

The one that matters most is the redundancy test. An early version of this code
produced correct, nearly complete contigs and a great deal of them, because at
thirty fold coverage every read is free to seed its own chain along the same
stretch of genome. A 120 kb genome came out as 2.8 Mb spread over a hundred
contigs, and every single contig was a perfectly good copy of the answer. Nothing
in the usual summary numbers makes that obvious, so it has its own test.
"""

from __future__ import annotations

import random

import pytest

from nagasm.olc import assemble_olc, overlap_summary
from nagasm.stats import summarise

from synthetic import random_sequence, reverse_complement, sample_reads


def genome_kmers(truth: str, k: int = 31) -> set[str]:
    forward = {truth[i : i + k] for i in range(len(truth) - k + 1)}
    return forward | {reverse_complement(x) for x in forward}


def accuracy(sequence: str, present: set[str], k: int = 31) -> float:
    if len(sequence) < k:
        return 0.0
    return sum(
        1 for i in range(len(sequence) - k + 1) if sequence[i : i + k] in present
    ) / (len(sequence) - k + 1)


def make_reads(length=120_000, depth=25, read_length=4_000, error_rate=0.0, seed=8):
    rng = random.Random(seed)
    truth = random_sequence(length, rng, gc=0.35)
    reads = sample_reads(
        truth, depth=depth, mean_length=read_length, rng=rng, error_rate=error_rate
    )
    return truth, reads


def test_a_clean_genome_comes_back_as_one_accurate_contig():
    truth, reads = make_reads()
    contigs = assemble_olc(reads)
    assert contigs, "something should be assembled"
    longest = max(contigs, key=lambda c: len(c.seq))
    assert len(longest.seq) > 0.9 * len(truth)
    assert accuracy(longest.seq, genome_kmers(truth)) == pytest.approx(1.0)


def test_the_assembly_is_not_the_genome_over_and_over():
    truth, reads = make_reads(depth=30)
    contigs = assemble_olc(reads)
    total = summarise(contigs).total_bp
    assert total < 1.5 * len(truth), (
        f"the assembly totals {total:,} bp for a {len(truth):,} bp genome, which "
        f"means reads are being used in more than one contig"
    )


def test_coverage_barely_changes_the_total_length():
    # Doubling the reads should give the same genome back, not twice as much of
    # it. This is the same redundancy problem seen from another angle.
    _, thin = make_reads(depth=15)
    _, thick = make_reads(depth=40)
    thin_total = summarise(assemble_olc(thin)).total_bp
    thick_total = summarise(assemble_olc(thick)).total_bp
    assert 0.7 < thick_total / thin_total < 1.4


def test_contigs_are_sorted_longest_first_and_named_in_order():
    _, reads = make_reads()
    contigs = assemble_olc(reads)
    lengths = [len(c.seq) for c in contigs]
    assert lengths == sorted(lengths, reverse=True)
    assert contigs[0].name == "contig_1"


def test_short_contigs_are_dropped():
    _, reads = make_reads()
    contigs = assemble_olc(reads, min_contig_bp=50_000)
    assert all(len(c.seq) >= 50_000 for c in contigs)


def test_no_reads_gives_no_contigs():
    assert assemble_olc([]) == []


def test_reads_that_overlap_nothing_still_come_out():
    rng = random.Random(4)
    unrelated = [random_sequence(3_000, rng) for _ in range(4)]
    contigs = assemble_olc(unrelated, min_contig_bp=500)
    assert len(contigs) == 4


def test_overlap_summary_reports_a_healthy_graph_for_deep_coverage():
    _, reads = make_reads(depth=25)
    info = overlap_summary(reads)
    assert info["n_reads"] == len(reads)
    assert info["share_isolated"] < 0.05
    assert info["mean_partners"] > 5


def test_overlap_summary_reports_isolation_when_reads_do_not_overlap():
    rng = random.Random(6)
    scattered = [random_sequence(2_000, rng) for _ in range(20)]
    info = overlap_summary(scattered)
    assert info["share_isolated"] > 0.9
