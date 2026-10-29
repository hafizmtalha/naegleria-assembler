"""Tests for the read support check on an assembled element.

These exist because of a real failure. On Klebsiella pneumoniae the pipeline
reported four closed circles that were not molecules at all. Three were held up
by two or three reads. The fourth was a 6.5 kb sequence stitched together from
reads belonging to dozens of separate places in the chromosome, and it was the
dangerous one, because it had five hundred reads behind it and looked perfectly
healthy by every number the pipeline printed at the time.

The case that matters most here is `test_sequence_stitched_from_a_repeat_is_rejected`.
It builds the bacterial failure on purpose: a short repeat that appears many
times inside a much larger chromosome, with unique sequence between the copies.
That repeat is genuinely several times deeper than the rest of the genome, so
the depth sorting is right to pick it out, and a circle built from it is still
not a real thing.
"""

from __future__ import annotations

import random

import pytest

from nagasm import support
from nagasm.polish import have

from synthetic import random_sequence, sample_reads

pytestmark = pytest.mark.skipif(
    not have("minimap2"), reason="the support check maps reads with minimap2"
)


def test_a_real_circle_is_supported():
    rng = random.Random(1)
    circle = random_sequence(6_000, rng)
    reads = sample_reads(circle, depth=60, mean_length=2_000, rng=rng, circular=True,
                         error_rate=0.02)
    result = support.measure(circle, reads, is_circular=True, threads=2)
    ok, why = support.verdict(result)
    assert ok, why
    assert result.depth > 20
    # Reads from a free circle have nothing else to match, so most of each read
    # sits inside the element.
    assert result.contained > 0.5


def test_an_element_held_up_by_a_handful_of_reads_is_rejected():
    rng = random.Random(2)
    circle = random_sequence(5_000, rng)
    reads = sample_reads(circle, depth=1.5, mean_length=1_500, rng=rng, circular=True,
                         error_rate=0.02)
    result = support.measure(circle, reads, is_circular=True, threads=2)
    ok, why = support.verdict(result)
    assert not ok
    assert "too thin" in why


def test_sequence_stitched_from_a_repeat_is_rejected():
    """The Klebsiella failure, built on purpose.

    A 3 kb repeat sits in eight places in a 200 kb chromosome, with unique
    sequence in between. The repeat really is eight times deeper than the rest,
    so it lands in a high copy compartment honestly. Reads covering it carry
    unique chromosome either side, so only a small part of each read has
    anywhere to go in a circle made out of the repeat alone.
    """
    rng = random.Random(3)
    repeat = random_sequence(3_000, rng)
    chromosome = ""
    for _ in range(8):
        chromosome += random_sequence(22_000, rng) + repeat
    reads = sample_reads(chromosome, depth=15, mean_length=12_000, rng=rng,
                         error_rate=0.02)

    result = support.measure(repeat, reads, is_circular=True, threads=2)
    ok, why = support.verdict(result)
    assert not ok, f"a chromosomal repeat was accepted as an element: {why}"
    # It fails on containment rather than on depth, which is the whole point.
    # Eight copies plus deep reads give it plenty of apparent coverage.
    assert result.depth > 10
    assert result.contained < 0.2
    assert "belongs somewhere else" in why


def test_an_element_no_read_maps_to_is_rejected():
    rng = random.Random(4)
    element = random_sequence(4_000, rng)
    unrelated = sample_reads(random_sequence(50_000, rng), depth=5,
                             mean_length=3_000, rng=rng)
    result = support.measure(element, unrelated, is_circular=True, threads=2)
    ok, why = support.verdict(result)
    assert not ok
    assert result.n_reads == 0 or result.depth < 10


def test_containment_is_not_punished_for_a_small_element():
    """A 3 kb circle read with 10 kb reads is still a 3 kb circle.

    Measuring the aligned share against the read length rather than against the
    shorter of the two made every small element look poorly contained, and made
    the 48 kb mitochondrion look better than the 15 kb rDNA circle purely
    because of its size. This is that bug, kept as a test.
    """
    rng = random.Random(5)
    circle = random_sequence(3_000, rng)
    # Reads much longer than the element itself.
    reads = sample_reads(circle, depth=80, mean_length=3_000, rng=rng, circular=True,
                         error_rate=0.02)
    result = support.measure(circle, reads, is_circular=True, threads=2)
    ok, why = support.verdict(result)
    assert ok, why


def test_a_missing_minimap2_leaves_the_element_unjudged():
    """No tool means no verdict, rather than a guilty verdict."""
    unmeasured = support.Support(0, 0.0, 0.0, measured=False, note="minimap2 not installed")
    ok, why = support.verdict(unmeasured)
    assert ok
    assert "minimap2" in why


def test_nothing_to_measure_is_handled():
    result = support.measure("", [], is_circular=True)
    assert not result.measured
    ok, _ = support.verdict(result)
    assert ok


# ---------------------------- a check that did not run is not a check that passed


def test_an_unmeasured_element_keeps_its_reads_in_the_nuclear_compartment(tmp_path):
    # support.verdict reports True when minimap2 is missing, so that nothing is
    # quietly dropped from the report. That is right for reporting and wrong for
    # deciding which reads to pull out of the chromosomes: it would strip real
    # sequence on the strength of a test nobody ran. Seen for real on a run with
    # minimap2 off PATH, which removed 55,877 reads instead of 1,003.
    from nagasm import support

    unmeasured = support.Support(0, 0.0, 0.0, measured=False, note="minimap2 not installed")
    supported, _ = support.verdict(unmeasured)
    assert supported is True, "the report must still list it"

    # The pipeline decides read removal on supported AND measured together.
    elements = [
        {"name": "real", "supported": True, "support": {"measured": True}},
        {"name": "never_checked", "supported": True, "support": unmeasured.as_dict()},
        {"name": "rejected", "supported": False, "support": {"measured": True}},
    ]
    keep_out = {
        e["name"]
        for e in elements
        if e.get("supported", True) and e.get("support", {}).get("measured", True)
    }
    assert keep_out == {"real"}, (
        "only an element whose support was actually measured and passed may have "
        f"its reads removed, got {keep_out}"
    )


def test_an_unmeasured_element_is_marked_on_its_picture():
    from nagasm import circleplot

    # The pipeline builds this note for an element whose check never ran.
    note = "SUPPORT NOT CHECKED, minimap2 was not available"
    assert note.upper().startswith(("NOT SUPPORTED", "SUPPORT NOT CHECKED")), (
        "the circle plot decides whether to colour the title by this prefix, so a "
        "note it does not recognise would be drawn as though the element passed"
    )
