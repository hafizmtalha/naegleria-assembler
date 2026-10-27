"""Tests for the optional circular plot module.

The pictures themselves cannot be checked by eye in a test, so these check the
numbers behind the rings, which is where a plot actually goes wrong. A plot
that draws a smooth ring of the wrong values looks fine and tells you a lie.
"""

from __future__ import annotations

import random

import numpy as np
import pytest

from nagasm import circleplot
from nagasm.polish import have
from nagasm.seqio import Record, write_fasta


def random_sequence(n: int, rng: random.Random, gc: float = 0.5) -> str:
    out = []
    for _ in range(n):
        if rng.random() < gc:
            out.append(rng.choice("GC"))
        else:
            out.append(rng.choice("AT"))
    return "".join(out)


def test_gc_content_follows_the_sequence():
    sequence = "G" * 500 + "A" * 500
    values = circleplot.gc_content(sequence, window=250)
    assert values[0] == pytest.approx(1.0)
    assert values[1] == pytest.approx(1.0)
    assert values[2] == pytest.approx(0.0)
    assert values[3] == pytest.approx(0.0)


def test_gc_skew_has_the_sign_the_strand_implies():
    # All G means every base is on the G rich strand, so the skew is +1.
    assert circleplot.gc_skew("G" * 200, window=200)[0] == pytest.approx(1.0)
    assert circleplot.gc_skew("C" * 200, window=200)[0] == pytest.approx(-1.0)
    # No G and no C at all is undefined, and must not divide by zero.
    assert circleplot.gc_skew("A" * 200, window=200)[0] == 0.0


def test_cumulative_skew_turns_where_the_strand_turns():
    # A sequence that is G rich for the first half and C rich for the second
    # is the shape a replicon makes. The running total should peak at the
    # switch, which is the middle.
    sequence = "G" * 1_000 + "C" * 1_000
    values = circleplot.cumulative_gc_skew(sequence, window=50)
    peak = int(np.argmax(values))
    assert 18 <= peak <= 21, f"the turn was found at window {peak}, not near 20"


def test_windows_cover_everything_exactly_once():
    bounds = circleplot._windows(1_050, 100)
    assert bounds[0] == (0, 100)
    assert bounds[-1] == (1_000, 1_050), "the short last window must not be dropped"
    assert sum(e - s for s, e in bounds) == 1_050


def test_a_plot_is_written_and_is_not_empty(tmp_path):
    rng = random.Random(1)
    sequence = random_sequence(6_000, rng)
    path = circleplot.plot_circle(sequence, tmp_path / "c.png", name="test_circle")
    assert path.exists()
    # A blank canvas comes out tiny. A real plot does not.
    assert path.stat().st_size > 20_000


def test_too_short_to_draw_is_refused(tmp_path):
    with pytest.raises(ValueError):
        circleplot.plot_circle("ACGT" * 10, tmp_path / "c.png")


def test_plot_all_draws_one_picture_per_record(tmp_path):
    rng = random.Random(2)
    fasta = tmp_path / "elements.fasta"
    write_fasta(
        fasta,
        [
            Record("first_circular", random_sequence(3_000, rng)),
            Record("second_linear", random_sequence(4_000, rng)),
            Record("far_too_short", "ACGT" * 5),
        ],
    )
    written = circleplot.plot_all(fasta, tmp_path / "out")
    assert len(written) == 2, "the short record should be skipped, not drawn"
    assert all(p.exists() for p in written)


def test_a_name_with_awkward_characters_is_still_a_valid_filename(tmp_path):
    rng = random.Random(5)
    fasta = tmp_path / "e.fasta"
    write_fasta(fasta, [Record("weird name/with:slashes", random_sequence(2_000, rng))])
    written = circleplot.plot_all(fasta, tmp_path / "out")
    assert len(written) == 1
    assert "/" not in written[0].name and ":" not in written[0].name


@pytest.mark.skipif(not have("minimap2"), reason="needs minimap2")
def test_the_depth_ring_finds_a_region_that_is_covered_twice(tmp_path):
    # A circle where one quarter is present in two copies in the read set
    # should come back with that quarter at roughly double depth.
    rng = random.Random(7)
    circle = random_sequence(4_000, rng)
    reads = []
    for _ in range(400):
        start = rng.randrange(0, 4_000)
        piece = (circle + circle)[start : start + 600]
        reads.append(piece)
    # Extra reads covering only the first quarter.
    for _ in range(400):
        start = rng.randrange(0, 400)
        reads.append(circle[start : start + 600])

    depth = circleplot.read_depth(circle, reads, window=200, threads=2)
    assert depth is not None
    first_quarter = depth[:4].mean()
    rest = depth[6:].mean()
    assert first_quarter > rest * 1.4, (
        f"the doubled region came back at {first_quarter:.1f} against {rest:.1f} "
        "elsewhere, so the depth ring is not tracking the reads"
    )


@pytest.mark.skipif(not have("minimap2"), reason="needs minimap2")
def test_reference_rings_report_high_identity_for_a_near_perfect_match(tmp_path):
    rng = random.Random(11)
    sequence = random_sequence(5_000, rng)
    reference = tmp_path / "ref.fasta"
    write_fasta(reference, [Record("ref", sequence)], wrap=0)

    rings = circleplot.reference_rings(sequence, reference, window=100)
    assert rings, "an exact match produced no reference rings"
    identity = rings[0]
    assert identity.kind == "heat"
    assert identity.values.mean() > 0.95
    # An exact match has nothing missing, so there is no second ring.
    assert len(rings) == 1


@pytest.mark.skipif(not have("minimap2"), reason="needs minimap2")
def test_sequence_absent_from_the_reference_is_reported_as_absent(tmp_path):
    rng = random.Random(13)
    shared = random_sequence(4_000, rng)
    unique = random_sequence(4_000, random.Random(99))
    reference = tmp_path / "ref.fasta"
    write_fasta(reference, [Record("ref", shared)], wrap=0)

    rings = circleplot.reference_rings(shared + unique, reference, window=100)
    assert len(rings) == 2, "the half that is not in the reference was not flagged"
    missing = rings[1]
    assert missing.label == "not in reference"
    # Roughly the second half should be flagged.
    assert 0.3 < missing.values.mean() < 0.7


def test_no_reference_match_is_an_answer_not_a_crash(tmp_path):
    if not have("minimap2"):
        pytest.skip("needs minimap2")
    rng = random.Random(17)
    reference = tmp_path / "ref.fasta"
    write_fasta(reference, [Record("ref", random_sequence(3_000, rng))], wrap=0)
    rings = circleplot.reference_rings(
        random_sequence(3_000, random.Random(18)), reference, window=100
    )
    assert rings == []


def test_depth_ring_without_minimap2_or_reads_is_skipped_quietly():
    assert circleplot.read_depth("ACGT" * 500, [], window=100) is None
