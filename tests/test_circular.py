"""Tests for the circular assembler.

Three of these exist because the code failed them while it was being written,
and each failure was invisible in the summary numbers. The assembly reported a
sensible length and claimed to be circular while the sequence itself was noise.
So the tests here check the sequence against a known answer, not just its size.

The failures worth knowing about were:

* the walk going round and round a circle forever, because coverage is deep
  enough that it meets fresh reads on every lap,
* reverse complemented reads being placed a k-mer's width out, because a shift
  relates the starts of two k-mers and flipping swaps which end that is,
* reads being attached to the layout using a shift written from the wrong read's
  point of view.
"""

from __future__ import annotations

import random

import pytest

from nagasm.circular import (
    _place_neighbour,
    assemble_circle,
    downsample,
    repeat_content,
    resolvability_note,
    rotate_to_canonical_start,
)

from synthetic import random_sequence, reverse_complement, sample_reads


def circular_identity(assembly: str, truth: str, k: int = 31) -> float:
    """Share of the assembly's k-mers found anywhere in the true circle.

    Both strands are tried, because a circle has no natural orientation and an
    assembly of one is equally right either way round. Checking only the forward
    strand makes a perfect assembly look like a total failure, which cost a while
    to spot the first time.
    """
    if not assembly:
        return 0.0
    best = 0.0
    for candidate in (truth, reverse_complement(truth)):
        doubled = candidate + candidate
        present = {doubled[i : i + k] for i in range(len(candidate))}
        hits = sum(1 for i in range(len(assembly) - k + 1) if assembly[i : i + k] in present)
        best = max(best, hits / max(len(assembly) - k + 1, 1))
    return best


def make_circle(length=12_000, depth=200, read_length=3_000, error_rate=0.0, seed=42):
    rng = random.Random(seed)
    truth = random_sequence(length, rng, gc=0.33)
    reads = sample_reads(
        truth, depth=depth, mean_length=read_length, rng=rng, circular=True,
        error_rate=error_rate,
    )
    return truth, reads


# --- placing one read against another ------------------------------------


@pytest.mark.parametrize("anchor_forward", [True, False])
@pytest.mark.parametrize("other_forward", [True, False])
def test_placement_is_exact_in_all_four_orientation_cases(anchor_forward, other_forward):
    # Two reads cut from one sequence at known places, each optionally flipped.
    # Whatever the combination, the placement has to put the second read back
    # where it came from.
    from nagasm.circular import _overlap_graph

    rng = random.Random(3)
    truth = random_sequence(20_000, rng, gc=0.35)
    k = 15
    a_start, b_start, span = 0, 1_500, 4_000

    a = truth[a_start : a_start + span]
    b = truth[b_start : b_start + span]
    reads = [
        a if anchor_forward else reverse_complement(a),
        b if other_forward else reverse_complement(b),
    ]

    graph = _overlap_graph(reads, k, 10, 1.0, 4, 3)
    assert graph.get(0), "the two reads overlap by 2.5 kb and must be found"

    other, same_strand, shift, _ = graph[0][0]
    assert other == 1
    offset, forward = _place_neighbour(
        a_start, anchor_forward, len(reads[0]), len(reads[1]), same_strand, shift, k
    )
    assert abs(offset - b_start) <= 2, f"placed at {offset}, should be {b_start}"
    # The anchor was handed to the placer already oriented so that the draft runs
    # the same way as the original sequence. So the second read has to be laid
    # down as it is exactly when it was cut out the right way round.
    assert forward == other_forward


# --- whole circles --------------------------------------------------------


def test_a_clean_circle_is_recovered_exactly():
    truth, reads = make_circle(error_rate=0.0)
    asm = assemble_circle(reads, expected_bp=len(truth), target_depth=50)

    assert asm.is_circular
    assert len(asm) == len(truth)
    assert circular_identity(asm.sequence, truth) == pytest.approx(1.0)


def test_the_walk_stops_after_one_lap():
    # Deep coverage means the walk meets fresh reads after going right round. If
    # it does not notice, it laps the circle again and again. This test fails
    # loudly if that ever comes back: the length runs into the hundreds of
    # kilobases for a twelve kilobase circle.
    truth, reads = make_circle(depth=400, error_rate=0.0)
    asm = assemble_circle(reads, expected_bp=len(truth), target_depth=120)
    assert len(asm) < len(truth) * 1.1, f"the walk lapped the circle: {len(asm)} bp"


def test_circularity_is_detected_and_the_length_is_right_despite_noise():
    for error_rate in (0.01, 0.02, 0.05):
        truth, reads = make_circle(error_rate=error_rate)
        asm = assemble_circle(reads, expected_bp=len(truth), target_depth=50)
        assert asm.is_circular, f"should close at {error_rate:.0%} error"
        assert abs(len(asm) - len(truth)) / len(truth) < 0.01


def test_a_noisy_draft_carries_read_error_and_says_so():
    # The draft is stitched raw read sequence, so its accuracy should be about
    # what the reads were. We check it is in the right ballpark rather than
    # pinning a number, and that the user is told it needs polishing.
    truth, reads = make_circle(error_rate=0.02)
    asm = assemble_circle(reads, expected_bp=len(truth), target_depth=50)
    identity = circular_identity(asm.sequence, truth)
    expected = 0.98**31  # a 31-mer survives only if every base in it is right
    assert 0.5 * expected < identity < 2.0 * expected
    assert any("polish" in note.lower() for note in asm.notes)


def test_a_linear_element_is_not_called_circular():
    rng = random.Random(9)
    truth = random_sequence(12_000, rng, gc=0.33)
    reads = sample_reads(truth, depth=80, mean_length=3_000, rng=rng, circular=False)
    asm = assemble_circle(reads, expected_bp=len(truth), target_depth=50)
    assert not asm.is_circular
    assert any("line" in note for note in asm.notes)


def test_too_few_reads_gives_an_empty_answer_rather_than_an_exception():
    asm = assemble_circle(["ACGT" * 20], expected_bp=1_000)
    assert asm.sequence == ""
    assert asm.notes


def test_no_reads_at_all_is_handled():
    asm = assemble_circle([], expected_bp=1_000)
    assert asm.sequence == ""
    assert not asm.is_circular


# --- the pieces around the assembler --------------------------------------


def test_rotation_makes_two_spins_of_one_circle_agree():
    rng = random.Random(11)
    circle = random_sequence(5_000, rng, gc=0.33)
    spun = circle[1_234:] + circle[:1_234]
    assert rotate_to_canonical_start(circle) == rotate_to_canonical_start(spun)


def test_rotation_keeps_the_sequence_the_same_length_and_content():
    rng = random.Random(12)
    circle = random_sequence(2_000, rng)
    rotated = rotate_to_canonical_start(circle)
    assert len(rotated) == len(circle)
    assert sorted(rotated) == sorted(circle)


def test_rotation_leaves_a_very_short_sequence_alone():
    assert rotate_to_canonical_start("ACGT", k=16) == "ACGT"


def test_downsampling_keeps_the_longest_reads_and_hits_the_budget():
    reads = [f"{'A' * n}" for n in (100, 5_000, 300, 4_000, 200)]
    kept = downsample(reads, target_depth=2, expected_bp=4_000)
    assert len(kept[0]) == 5_000
    assert sum(len(r) for r in kept) >= 8_000


def test_downsampling_with_no_expected_size_keeps_everything():
    reads = ["ACGT" * 100, "TTTT" * 100]
    assert len(downsample(reads, 10, 0)) == 2


def test_repeat_content_spots_a_tandem_repeat_and_ignores_unique_sequence():
    rng = random.Random(13)
    unique = random_sequence(4_000, rng)
    assert repeat_content(unique, k=50) < 0.05

    block = random_sequence(500, rng)
    repetitive = (block * 6) + unique[:1_000]
    assert repeat_content(repetitive, k=50) > 0.5


def test_the_resolvability_note_warns_when_reads_are_shorter_than_the_repeats():
    rng = random.Random(14)
    block = random_sequence(900, rng)
    element = block * 8
    short_reads = ["A" * 400]
    note = resolvability_note(element, short_reads)
    assert "unsettled" in note

    long_reads = ["A" * (len(element) + 100)]
    assert "resolved" in resolvability_note(element, long_reads)
