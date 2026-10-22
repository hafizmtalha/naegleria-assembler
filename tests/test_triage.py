"""Tests for the compartment triage.

The point of most of these is not that some number comes out right. It is that
the triage keeps working on a read set shaped like Naegleria, where the
compartment we care about most is a thousandth of the size of the biggest one.
Several obvious ways of writing this code pass on an ordinary genome and fail
here, so there is a test for each of those ways.
"""

from __future__ import annotations

import collections

import numpy as np
import pytest

from nagasm.triage import (
    Compartment,
    find_compartments,
    find_error_threshold,
    triage,
)

from synthetic import naegleria_like, random_sequence, sample_reads


def purity(result, labels) -> float:
    """What share of placed reads sit with the majority of their own kind."""
    correct = total = 0
    for i in range(len(result.compartments)):
        counts = collections.Counter(
            lab for lab, sel in zip(labels, result.assignment == i) if sel
        )
        if counts:
            correct += max(counts.values())
            total += sum(counts.values())
    return correct / total if total else 0.0


def compartment_for(result, labels, source: str) -> Compartment | None:
    """The compartment holding most reads from a given true source."""
    counts = collections.Counter(
        int(a) for a, lab in zip(result.assignment, labels) if lab == source and a >= 0
    )
    if not counts:
        return None
    return result.compartments[counts.most_common(1)[0][0]]


# --- the error threshold -------------------------------------------------


def test_error_threshold_sits_in_the_dip_after_the_error_spike():
    # A steep fall from a count of one, a dip at 4, then a real hump.
    hist = np.array([0, 900, 300, 90, 20, 60, 160, 300, 160, 60, 20, 5], dtype=np.int64)
    assert find_error_threshold(hist) == 4


def test_error_threshold_gives_up_rather_than_inventing_a_cut():
    # A histogram that only ever falls has no dip, which is what low coverage
    # looks like. Returning 1 keeps every k-mer instead of guessing.
    hist = np.array([0, 900, 300, 90, 30, 10, 3, 1], dtype=np.int64)
    assert find_error_threshold(hist) == 1


def test_error_threshold_copes_with_a_tiny_histogram():
    assert find_error_threshold(np.array([0, 5], dtype=np.int64)) == 1


# --- picking compartments out of a histogram -----------------------------


def _histogram_with_humps(humps, size=20_001):
    """Build a k-mer histogram with Gaussian humps at given depth and area."""
    hist = np.zeros(size, dtype=float)
    depths = np.arange(size)
    for peak, area in humps:
        width = max(peak * 0.18, 1.0)
        hist += area * np.exp(-0.5 * ((depths - peak) / width) ** 2)
    return hist.astype(np.int64)


def test_a_tiny_deep_compartment_is_not_lost_next_to_a_huge_shallow_one():
    # This is the Naegleria case. The second hump covers about a thousandth of
    # the sequence of the first, so its peak is far shorter. Judging peaks by
    # height relative to the tallest would throw it away, which is the single
    # easiest way to get this organism wrong.
    hist = _histogram_with_humps([(60, 40_000), (4_000, 40)])
    comps, nuclear = find_compartments(hist, scale=1)

    assert len(comps) == 2, f"both humps should survive, got {[c.label for c in comps]}"
    assert comps[nuclear].label == "nuclear"
    assert comps[nuclear].peak_depth < 200
    deep = [c for c in comps if c.label != "nuclear"]
    assert deep and deep[0].peak_depth > 1_000


def test_the_nuclear_compartment_is_the_widest_not_the_tallest():
    # Here the deep hump is given a larger area than the shallow one, as happens
    # when a high copy element contributes more reads than the chromosomes.
    # Nuclear must still be chosen by span, which the depth axis weighting makes
    # possible.
    hist = _histogram_with_humps([(50, 30_000), (3_000, 300)])
    comps, nuclear = find_compartments(hist, scale=1)
    assert comps[nuclear].peak_depth < 200
    assert comps[nuclear].span_bp == max(c.span_bp for c in comps)


def test_span_estimates_land_near_the_truth():
    hist = _histogram_with_humps([(60, 40_000), (4_000, 40)])
    comps, nuclear = find_compartments(hist, scale=1)
    # The area under a Gaussian is amplitude times width times root two pi.
    expected_nuclear = 40_000 * (60 * 0.18) * np.sqrt(2 * np.pi)
    got = comps[nuclear].span_bp
    assert 0.5 * expected_nuclear < got < 2.0 * expected_nuclear


def test_scale_multiplies_through_to_the_span():
    hist = _histogram_with_humps([(60, 40_000), (4_000, 40)])
    one, _ = find_compartments(hist, scale=1)
    hundred, _ = find_compartments(hist, scale=100)
    assert len(one) == len(hundred)
    for a, b in zip(one, hundred):
        assert b.span_bp == a.span_bp * 100


def test_a_featureless_histogram_still_gives_one_compartment():
    hist = np.zeros(2_001, dtype=np.int64)
    hist[40:60] = 1_000
    comps, nuclear = find_compartments(hist, scale=1)
    assert len(comps) >= 1
    assert nuclear >= 0


def test_an_empty_histogram_gives_nothing_rather_than_raising():
    comps, nuclear = find_compartments(np.zeros(5, dtype=np.int64), scale=1)
    assert comps == [] or nuclear >= 0


# --- end to end on synthetic reads ---------------------------------------


def test_clean_reads_separate_perfectly():
    g = naegleria_like(seed=1, error_rate=0.0, nuclear_depth=60, mito_depth=700, rdna_depth=4000)
    res = triage(g.reads, k=21, scale=20)

    assert purity(res, g.labels) > 0.99
    assert res.nuclear is not None
    assert len(res.compartments) == 3

    nuc = compartment_for(res, g.labels, "nuclear")
    mito = compartment_for(res, g.labels, "mito")
    rdna = compartment_for(res, g.labels, "rdna")
    assert nuc.label == "nuclear"
    assert mito.peak_depth > nuc.peak_depth
    assert rdna.peak_depth > mito.peak_depth


def test_the_high_copy_elements_are_reported_as_small():
    g = naegleria_like(seed=1, error_rate=0.0, nuclear_depth=60, mito_depth=700, rdna_depth=4000)
    res = triage(g.reads, k=21, scale=20)
    for c in res.high_copy:
        assert c.span_bp < res.nuclear.span_bp / 10


def test_a_deep_compartment_with_more_reads_does_not_steal_the_nuclear_label():
    # The rDNA circle here contributes more total sequence than the chromosomes,
    # which is the situation that a read count based rule gets wrong.
    g = naegleria_like(
        seed=5,
        nuclear_bp=200_000,
        nuclear_depth=20,
        rdna_bp=4_000,
        rdna_depth=3_000,
        mito_bp=12_000,
        mito_depth=200,
    )
    by_source = collections.Counter()
    for read, lab in zip(g.reads, g.labels):
        by_source[lab] += len(read)
    assert by_source["rdna"] > by_source["nuclear"], "the test setup should be lopsided"

    res = triage(g.reads, k=21, scale=20)
    nuc = compartment_for(res, g.labels, "nuclear")
    assert nuc is not None and nuc.label == "nuclear"


@pytest.mark.parametrize("error_rate", [0.0, 0.02, 0.05])
def test_separation_survives_realistic_long_read_error(error_rate):
    # Modern nanopore and HiFi sit in this range. Above it the method degrades,
    # which is covered by its own test below.
    g = naegleria_like(
        seed=4, error_rate=error_rate, nuclear_depth=60, mito_depth=700, rdna_depth=4000
    )
    res = triage(g.reads, k=21, scale=20)
    assert purity(res, g.labels) > 0.95
    assert (res.assignment >= 0).mean() > 0.95


def test_very_noisy_reads_degrade_with_a_warning_rather_than_silently():
    g = naegleria_like(
        seed=4, error_rate=0.20, nuclear_depth=60, mito_depth=700, rdna_depth=4000
    )
    res = triage(g.reads, k=21, scale=20)
    # We do not insist the answer is good at this error rate. We insist the
    # user is told, because a quiet wrong answer is the worst outcome.
    assert res.warnings, "a read set this noisy should produce warnings"


def test_homopolymer_compression_is_chosen_for_long_reads_and_not_for_short():
    import random

    long_reads = naegleria_like(seed=2, mean_read_length=3000).reads[:600]
    assert triage(long_reads, k=21, scale=20).hpc is True

    rng = random.Random(0)
    short_reads = [random_sequence(150, rng) for _ in range(800)]
    assert triage(short_reads, k=15, scale=5).hpc is False


def test_read_indices_round_trips_the_labels():
    g = naegleria_like(seed=1, error_rate=0.0)
    res = triage(g.reads, k=21, scale=20)
    total = 0
    for c in res.compartments:
        idx = res.read_indices(c.label)
        assert len(idx) == c.n_reads
        total += len(idx)
    assert total == int((res.assignment >= 0).sum())


def test_asking_for_a_compartment_that_is_not_there_raises():
    g = naegleria_like(seed=1, error_rate=0.0)
    res = triage(g.reads, k=21, scale=20)
    with pytest.raises(KeyError):
        res.read_indices("no_such_compartment")


def test_no_reads_is_an_error_not_a_crash():
    with pytest.raises(ValueError, match="no reads"):
        triage([])


def test_a_single_compartment_genome_is_handled():
    # Not every organism is Naegleria. A plain genome with one compartment must
    # not be split into imaginary ones.
    import random

    rng = random.Random(7)
    template = random_sequence(150_000, rng, gc=0.4)
    reads = sample_reads(template, depth=40, mean_length=2000, rng=rng)
    res = triage(reads, k=21, scale=20)
    assert res.nuclear is not None
    assert len(res.high_copy) == 0, "a plain genome should have no high copy compartments"


def test_total_bases_and_read_counts_add_up():
    g = naegleria_like(seed=3, error_rate=0.0)
    res = triage(g.reads, k=21, scale=20)
    assert sum(c.n_reads for c in res.compartments) == int((res.assignment >= 0).sum())
    placed = res.assignment >= 0
    assert sum(c.total_bases for c in res.compartments) == int(res.read_length[placed].sum())
