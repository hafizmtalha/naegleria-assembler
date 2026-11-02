"""Tests for the script that builds the simulated read set.

This script is how anyone running the examples gets data, so if it quietly
produces something the pipeline cannot use then nothing in the README works and
the failure looks like the pipeline's fault. It had two bugs that behaved exactly
that way, both of which exited cleanly, and both are pinned here.

The tests do not assemble anything. They check the shape of what comes out, since
that is what decides whether a run can demonstrate anything.
"""

from __future__ import annotations

import gzip
import importlib.util
import random
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]


def load_script():
    """Import the script by path, since scripts/ is not a package."""
    path = ROOT / "scripts" / "make_test_data.py"
    spec = importlib.util.spec_from_file_location("make_test_data", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


mtd = load_script()


def run(tmp_path, *args) -> int:
    """Run the script's main with argv set, the way a person would call it."""
    argv = ["make_test_data.py", "--out", str(tmp_path), *args]
    old = sys.argv
    sys.argv = argv
    try:
        return mtd.main()
    finally:
        sys.argv = old


def read_fastq(path: Path) -> list[str]:
    with gzip.open(path, "rt") as handle:
        lines = handle.read().splitlines()
    return lines[1::4]


# ------------------------------------------------------------- the pieces


def test_sequence_lands_near_the_gc_asked_for():
    rng = random.Random(1)
    seq = mtd.random_sequence(20_000, rng, gc=0.36)
    gc = sum(1 for b in seq if b in "GC") / len(seq)
    assert 0.34 < gc < 0.38, f"asked for 0.36 GC, got {gc:.3f}"


def test_reverse_complement_is_its_own_inverse():
    rng = random.Random(2)
    seq = mtd.random_sequence(500, rng)
    assert mtd.reverse_complement(mtd.reverse_complement(seq)) == seq
    assert mtd.reverse_complement("ACGT") == "ACGT"
    assert mtd.reverse_complement("AAAA") == "TTTT"


def test_repeats_go_in_the_right_number_of_times():
    rng = random.Random(3)
    backbone = mtd.random_sequence(50_000, rng)
    repeat = mtd.random_sequence(1_200, rng)
    out, at = mtd.with_repeats(backbone, repeat, 8, rng)
    assert len(at) == 8
    assert out.count(repeat) == 8, "a copy went missing, or two landed on top of each other"
    assert len(out) == len(backbone) + 8 * len(repeat)
    assert at == sorted(at), "the copies should be in order along the sequence"


def test_no_base_of_the_backbone_is_lost_to_the_repeats():
    rng = random.Random(4)
    backbone = mtd.random_sequence(30_000, rng)
    repeat = mtd.random_sequence(900, rng)
    out, _ = mtd.with_repeats(backbone, repeat, 5, rng)
    assert out.replace(repeat, "") == backbone


def test_errors_land_near_the_rate_asked_for():
    rng = random.Random(5)
    seq = mtd.random_sequence(40_000, rng)
    clean = mtd.mutate(seq, 0.0, rng)
    assert clean == seq, "a rate of zero must change nothing"
    noisy = mtd.mutate(seq, 0.10, rng)
    # Indels move the length, so compare lengths rather than position by position.
    drift = abs(len(noisy) - len(seq)) / len(seq)
    assert drift < 0.05, f"length moved by {drift:.1%}, which is more than indels explain"
    assert noisy != seq


def test_reads_from_a_circle_cross_the_join():
    # Without doubling the sequence first, every read stops at the same base, the
    # join is never covered, and the circle finder has nothing to close.
    rng = random.Random(6)
    circle = mtd.random_sequence(4_000, rng)
    reads = mtd.sample_reads(circle, 60, 1_500, 0.0, rng, circular=True)
    crossing = [r for r in reads if r not in circle and mtd.reverse_complement(r) not in circle]
    assert crossing, "no read ran across the join, so a circle could not be closed from these"


def test_reads_stay_inside_a_linear_molecule():
    rng = random.Random(7)
    seq = mtd.random_sequence(20_000, rng)
    reads = mtd.sample_reads(seq, 10, 2_000, 0.0, rng, circular=False)
    for read in reads:
        assert read in seq or mtd.reverse_complement(read) in seq


# ------------------------------------------------------- the whole script


def test_a_default_run_writes_everything_it_promises(tmp_path):
    assert run(tmp_path) == 0
    for name in ("reads.fastq.gz", "truth_chromosome.fasta", "truth_element.fasta",
                 "truth_repeat.fasta", "README.md"):
        assert (tmp_path / name).exists(), f"{name} was not written"
    reads = read_fastq(tmp_path / "reads.fastq.gz")
    assert len(reads) > 1_000
    assert all(len(r) >= 200 for r in reads), "a read came out too short to be useful"


def test_the_gzip_header_carries_no_timestamp(tmp_path):
    """The exact version of the test below, which on its own is unreliable.

    gzip puts the current time in bytes 4 to 7 of its header. Comparing two runs
    only notices that if the two happen to straddle a second boundary, and each
    run here takes a fraction of a second, so most of the time the stamps match
    and the comparison passes with the bug still in place. Reading the field
    directly has no such luck in it.
    """
    run(tmp_path, "--small")
    header = (tmp_path / "reads.fastq.gz").read_bytes()[:8]
    stamp = int.from_bytes(header[4:8], "little")
    assert stamp == 0, (
        f"the gzip header carries a timestamp of {stamp}, so the same seed will not "
        "give the same bytes twice"
    )


def test_the_same_seed_gives_the_same_file_byte_for_byte(tmp_path):
    # The promise the script makes to anyone reading its output. Kept because it
    # is what a person would check, with the exact version above to back it up.
    a, b = tmp_path / "a", tmp_path / "b"
    run(a, "--small")
    run(b, "--small")
    assert (a / "reads.fastq.gz").read_bytes() == (b / "reads.fastq.gz").read_bytes()
    assert (a / "truth_element.fasta").read_text() == (b / "truth_element.fasta").read_text()


def test_a_different_seed_gives_different_reads(tmp_path):
    a, b = tmp_path / "a", tmp_path / "b"
    run(a, "--small", "--seed", "1")
    run(b, "--small", "--seed", "2")
    assert (a / "reads.fastq.gz").read_bytes() != (b / "reads.fastq.gz").read_bytes()


def test_the_element_is_not_a_piece_of_the_chromosome(tmp_path):
    # If the circle were cut out of the chromosome, finding it would prove
    # nothing, and a match to it later would be meaningless.
    run(tmp_path, "--small")
    chromosome = "".join(
        l.strip() for l in (tmp_path / "truth_chromosome.fasta").read_text().splitlines()
        if not l.startswith(">")
    )
    element = "".join(
        l.strip() for l in (tmp_path / "truth_element.fasta").read_text().splitlines()
        if not l.startswith(">")
    )
    assert element not in chromosome
    assert element[:200] not in chromosome


def test_the_repeat_really_is_in_the_chromosome_several_times(tmp_path):
    run(tmp_path, "--small")
    chromosome = "".join(
        l.strip() for l in (tmp_path / "truth_chromosome.fasta").read_text().splitlines()
        if not l.startswith(">")
    )
    repeat = "".join(
        l.strip() for l in (tmp_path / "truth_repeat.fasta").read_text().splitlines()
        if not l.startswith(">")
    )
    assert chromosome.count(repeat) == 8


# ------------------------------------------- the two settings that broke


def test_a_flat_depth_spectrum_is_refused(tmp_path):
    # With the circle at the same depth as the chromosome there is no hump to
    # find, so the data would demonstrate nothing. Better to refuse than to hand
    # someone a set that cannot work.
    with pytest.raises(SystemExit) as exc:
        run(tmp_path, "--element-depth", "50", "--depth", "35")
    assert exc.value.code != 0


def test_small_keeps_the_chromosome_as_the_bigger_read_population(tmp_path):
    # This is the bug --small had. Shrinking the chromosome while leaving the
    # circle at 600x made the circle supply most of the reads, so the depth
    # sorting took the circle for the nuclear peak and found no element at all.
    # The run still exited cleanly, which is what made it worth a test.
    run(tmp_path, "--small")
    reads = read_fastq(tmp_path / "reads.fastq.gz")
    chromosome = "".join(
        l.strip() for l in (tmp_path / "truth_chromosome.fasta").read_text().splitlines()
        if not l.startswith(">")
    )
    element = "".join(
        l.strip() for l in (tmp_path / "truth_element.fasta").read_text().splitlines()
        if not l.startswith(">")
    )
    # Reads are noisy, so place each one by which source its first clean stretch
    # matches. A 40 base window is long enough to be unique to one of the two.
    from_chromosome = from_element = 0
    for read in reads:
        probe = read[:40]
        rc = mtd.reverse_complement(read)[:40]
        if probe in chromosome or rc in chromosome:
            from_chromosome += 1
        elif probe in element + element or rc in element + element:
            from_element += 1
    assert from_chromosome > from_element, (
        f"{from_element} reads came off the circle against {from_chromosome} off the "
        "chromosome. The chromosome has to be the bigger population or the depth "
        "sorting calls the circle nuclear and no element is ever found"
    )


def test_small_still_leaves_the_circle_clearly_deeper(tmp_path):
    # The other half of the same balance. The circle has to stand well above the
    # chromosome or there is no hump, and the script's own guard requires 4x.
    run(tmp_path, "--small")
    text = (tmp_path / "README.md").read_text()
    assert "30x" in text and "200x" in text
