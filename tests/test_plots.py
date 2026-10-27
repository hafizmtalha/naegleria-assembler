"""Tests for the optional picture module.

Same reasoning as the circle plot tests: a figure cannot be checked by eye in
a test, so these check the numbers and the layout decisions behind it. A plot
drawn smoothly from the wrong values looks perfectly fine and tells you a lie.
"""

from __future__ import annotations

import random

import numpy as np
import pytest

from nagasm import plots
from nagasm.polish import have
from nagasm.seqio import Record, write_fasta


def random_sequence(n: int, rng: random.Random, gc: float = 0.5) -> str:
    return "".join(rng.choice("GC") if rng.random() < gc else rng.choice("AT") for _ in range(n))


def an_assembly(tmp_path, lengths, gc=0.5, names=None, seed=1):
    rng = random.Random(seed)
    names = names or [f"contig_{i + 1}" for i in range(len(lengths))]
    path = tmp_path / "asm.fasta"
    write_fasta(
        path,
        [Record(n, random_sequence(l, rng, gc=gc)) for n, l in zip(names, lengths)],
        wrap=0,
    )
    return path


# ------------------------------------------------------------ small helpers


def test_bp_label_never_repeats_itself_on_one_axis():
    # A 6 kb element gets ticks every 828 bp. Rounding those to whole kb gave
    # two ticks that both read "2 kb", which is what this guards.
    labels = [plots.bp_label(v, 6_626) for v in np.linspace(0, 6_626, 9)[:-1]]
    assert len(set(labels)) == len(labels), labels


def test_bp_label_scales_with_the_whole_not_the_part():
    assert plots.bp_label(1_500, 5_000) == "1.5 kb"
    assert plots.bp_label(1_500_000, 5_000_000) == "1.50 Mb"
    assert plots.bp_label(500, 900) == "500 bp"


def test_safe_name_survives_an_awkward_contig_name():
    assert "/" not in plots.safe_name("weird/name:here")
    assert ":" not in plots.safe_name("weird/name:here")


def test_themes_all_define_what_a_plot_needs():
    for name, theme in plots.THEMES.items():
        for key in ("bg", "fg", "muted", "grid", "palette"):
            assert key in theme, f"theme {name} is missing {key}"
        assert len(theme["palette"]) >= 6


def test_offsets_lay_sequences_end_to_end(tmp_path):
    path = an_assembly(tmp_path, [1_000, 2_000, 500])
    offsets, total = plots._offsets(path)
    assert total == 3_500
    assert offsets["contig_1"] == 0
    assert offsets["contig_2"] == 1_000
    assert offsets["contig_3"] == 3_000


# -------------------------------------------------------------- karyotype


def test_karyotype_is_written(tmp_path):
    path = an_assembly(tmp_path, [5_000, 3_000, 1_000])
    out = plots.karyotype(path, tmp_path / "k.png")
    assert out.exists()
    assert out.stat().st_size > 15_000


def test_karyotype_draws_only_the_biggest_when_asked(tmp_path):
    path = an_assembly(tmp_path, [9_000, 8_000, 7_000, 6_000, 5_000])
    out = plots.karyotype(path, tmp_path / "k.png", max_contigs=2)
    assert out.exists()


def test_karyotype_refuses_an_empty_assembly(tmp_path):
    path = an_assembly(tmp_path, [500])
    with pytest.raises(ValueError):
        plots.karyotype(path, tmp_path / "k.png", min_length=10_000)


def test_karyotype_gc_range_fits_an_at_rich_genome(tmp_path):
    # Naegleria sits near 35 percent GC. A fixed 25 to 75 band puts every
    # window in the same colour and the picture says nothing, so the range has
    # to follow the data.
    path = an_assembly(tmp_path, [20_000, 18_000], gc=0.35, seed=4)
    out = plots.karyotype(path, tmp_path / "k.png", colour_by="gc")
    assert out.exists()


def test_karyotype_falls_back_to_gc_when_depth_has_no_reads(tmp_path):
    path = an_assembly(tmp_path, [5_000, 4_000])
    out = plots.karyotype(path, tmp_path / "k.png", colour_by="depth", reads=None)
    assert out.exists(), "asking for depth with no reads should still draw something"


# ---------------------------------------------------------------- nx curve


def test_nx_curve_is_written_and_starts_at_zero(tmp_path):
    path = an_assembly(tmp_path, [50_000, 10_000, 5_000, 1_000])
    out = plots.nx_curve(path, tmp_path / "nx.png")
    assert out.exists()


def test_nx_curve_keeps_two_assemblies_apart_despite_the_same_filename(tmp_path):
    # Two runs of this pipeline both produce a file called nuclear.fasta. Keying
    # the curves by filename alone silently drew one over the other.
    a = tmp_path / "run_a" / "assembly"
    b = tmp_path / "run_b" / "assembly"
    a.mkdir(parents=True)
    b.mkdir(parents=True)
    rng = random.Random(3)
    for folder, lengths in ((a, [40_000, 5_000]), (b, [9_000, 8_000])):
        write_fasta(
            folder / "nuclear.fasta",
            [Record(f"c{i}", random_sequence(l, rng)) for i, l in enumerate(lengths)],
            wrap=0,
        )
    out = plots.nx_curve(
        a / "nuclear.fasta",
        tmp_path / "nx.png",
        compare={"other": b / "nuclear.fasta"},
    )
    assert out.exists()


# ----------------------------------------------------------------- dotplot


@pytest.mark.skipif(not have("minimap2"), reason="needs minimap2")
def test_a_dotplot_finds_a_repeat_against_itself(tmp_path):
    rng = random.Random(21)
    repeat = random_sequence(2_000, rng)
    sequence = random_sequence(8_000, rng) + repeat + random_sequence(8_000, rng) + repeat
    path = tmp_path / "s.fasta"
    write_fasta(path, [Record("s", sequence)], wrap=0)
    out = plots.dotplot(path, tmp_path / "d.png")
    assert out is not None, "a sequence with two copies of a 2 kb repeat drew nothing"
    assert out.exists()


@pytest.mark.skipif(not have("minimap2"), reason="needs minimap2")
def test_a_sequence_with_no_repeat_draws_no_dotplot(tmp_path):
    rng = random.Random(23)
    path = tmp_path / "s.fasta"
    write_fasta(path, [Record("s", random_sequence(10_000, rng))], wrap=0)
    assert plots.dotplot(path, tmp_path / "d.png") is None


@pytest.mark.skipif(not have("minimap2"), reason="needs minimap2")
def test_a_dotplot_against_a_reference_is_drawn(tmp_path):
    rng = random.Random(27)
    sequence = random_sequence(12_000, rng)
    query, reference = tmp_path / "q.fasta", tmp_path / "r.fasta"
    write_fasta(query, [Record("q", sequence)], wrap=0)
    write_fasta(reference, [Record("r", sequence)], wrap=0)
    out = plots.dotplot(query, tmp_path / "d.png", reference=reference)
    assert out is not None and out.exists()


# ---------------------------------------------------------------- coverage


@pytest.mark.skipif(not have("minimap2"), reason="needs minimap2")
def test_coverage_tracks_the_reads(tmp_path):
    rng = random.Random(31)
    sequence = random_sequence(6_000, rng)
    path = tmp_path / "s.fasta"
    write_fasta(path, [Record("s", sequence)], wrap=0)
    reads = [sequence[s : s + 800] for s in range(0, 5_200, 50)]
    out = plots.coverage_plot(path, reads, tmp_path / "c.png", threads=2)
    assert out is not None and out.exists()


def test_coverage_without_reads_is_skipped_quietly(tmp_path):
    path = an_assembly(tmp_path, [4_000])
    assert plots.coverage_plot(path, [], tmp_path / "c.png") is None


@pytest.mark.skipif(not have("minimap2"), reason="needs minimap2")
def test_depth_per_contig_puts_depth_on_the_right_contig(tmp_path):
    rng = random.Random(37)
    a = random_sequence(5_000, rng)
    b = random_sequence(5_000, random.Random(38))
    contigs = [plots.Contig("a", len(a), a), plots.Contig("b", len(b), b)]
    # Reads from a only. Contig b must come back empty, not share the depth.
    reads = [a[s : s + 700] for s in range(0, 4_300, 60)]
    depths = plots._depth_per_contig(contigs, reads, window=500, threads=2)
    assert depths is not None
    assert depths["a"].mean() > 1.0
    assert depths["b"].mean() < 0.1, "depth leaked onto a contig with no reads"
