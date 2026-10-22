"""Tests for reading sequence files and for the assembly statistics."""

from __future__ import annotations

import gzip
import random

import pytest

from nagasm.seqio import (
    Record,
    read_fasta,
    read_fastq,
    read_sequences,
    reverse_complement,
    write_fasta,
    write_fastq,
)
from nagasm.stats import (
    coverage_estimate,
    gc_and_n,
    nx,
    read_length_summary,
    summarise,
)

from synthetic import random_sequence


# --- reading and writing --------------------------------------------------


def test_fasta_round_trips(tmp_path):
    records = [Record("one", "ACGTACGT"), Record("two", "TTTTGGGGCCCC")]
    path = tmp_path / "x.fasta"
    assert write_fasta(path, records) == 2
    back = list(read_fasta(path))
    assert [r.name for r in back] == ["one", "two"]
    assert [r.seq for r in back] == ["ACGTACGT", "TTTTGGGGCCCC"]


def test_wrapped_fasta_reads_the_same_as_unwrapped(tmp_path):
    seq = random_sequence(500, random.Random(1))
    wrapped = tmp_path / "w.fasta"
    flat = tmp_path / "f.fasta"
    write_fasta(wrapped, [Record("s", seq)], wrap=60)
    write_fasta(flat, [Record("s", seq)], wrap=0)
    assert list(read_fasta(wrapped))[0].seq == list(read_fasta(flat))[0].seq == seq


def test_fasta_names_stop_at_the_first_space(tmp_path):
    path = tmp_path / "d.fasta"
    path.write_text(">contig_1 length=500 coverage=30\nACGT\n")
    assert list(read_fasta(path))[0].name == "contig_1"


def test_fastq_round_trips_with_quality(tmp_path):
    records = [Record("r1", "ACGT", "IIII"), Record("r2", "TTGG", "####")]
    path = tmp_path / "x.fastq"
    write_fastq(path, records)
    back = list(read_fastq(path))
    assert [r.seq for r in back] == ["ACGT", "TTGG"]
    assert [r.qual for r in back] == ["IIII", "####"]


def test_fastq_without_quality_gets_a_placeholder(tmp_path):
    path = tmp_path / "x.fastq"
    write_fastq(path, [Record("r1", "ACGTA")])
    back = list(read_fastq(path))[0]
    assert len(back.qual) == len(back.seq)


def test_a_fastq_with_mismatched_quality_length_is_rejected(tmp_path):
    path = tmp_path / "bad.fastq"
    path.write_text("@r1\nACGTACGT\n+\nIII\n")
    with pytest.raises(ValueError, match="quality values"):
        list(read_fastq(path))


def test_gzipped_files_are_read_without_being_asked(tmp_path):
    path = tmp_path / "x.fasta.gz"
    with gzip.open(path, "wt") as handle:
        handle.write(">r1\nACGTACGT\n")
    assert list(read_sequences(path))[0].seq == "ACGTACGT"


def test_read_sequences_picks_the_format_by_itself(tmp_path):
    fa = tmp_path / "a.fa"
    fa.write_text(">r\nACGT\n")
    fq = tmp_path / "a.fq"
    fq.write_text("@r\nACGT\n+\nIIII\n")
    assert list(read_sequences(fa))[0].seq == "ACGT"
    assert list(read_sequences(fq))[0].seq == "ACGT"


def test_a_file_that_is_neither_format_is_rejected(tmp_path):
    path = tmp_path / "notes.txt"
    path.write_text("this is not a sequence file\n")
    with pytest.raises(ValueError, match="neither FASTA nor FASTQ"):
        list(read_sequences(path))


def test_an_empty_file_gives_no_records(tmp_path):
    path = tmp_path / "empty.fasta"
    path.write_text("")
    assert list(read_sequences(path)) == []


def test_reverse_complement_is_its_own_inverse():
    seq = random_sequence(300, random.Random(2))
    assert reverse_complement(reverse_complement(seq)) == seq


def test_reverse_complement_leaves_unknown_bases_alone():
    assert reverse_complement("ACGTN") == "NACGT"


# --- statistics -----------------------------------------------------------


def test_n50_is_the_length_that_half_the_assembly_sits_above():
    # Total 100, so the target is 50. Sorted longest first that is 50, 30, 15, 5.
    # The first contig alone reaches the target, so N50 is 50 and L50 is 1.
    assert nx([5, 50, 15, 30], 0.5) == (50, 1)

    # Total 100 again, but now no single contig reaches half. The running total
    # passes 50 at the second one, so N50 is 30 and L50 is 2.
    assert nx([5, 40, 25, 30], 0.5) == (30, 2)


def test_n90_reaches_further_down_the_list_than_n50():
    lengths = [100, 50, 20, 10, 5]
    n50, l50 = nx(lengths, 0.5)
    n90, l90 = nx(lengths, 0.9)
    assert n90 <= n50
    assert l90 >= l50


def test_nx_handles_an_empty_assembly():
    assert nx([], 0.5) == (0, 0)


def test_one_contig_is_its_own_n50():
    assert nx([1000], 0.5) == (1000, 1)


def test_gc_ignores_unknown_bases_but_counts_them_separately():
    gc, n = gc_and_n(["GGCCAATT", "NNNN"])
    assert gc == pytest.approx(50.0)
    assert n == pytest.approx(100.0 * 4 / 12)


def test_summarise_describes_an_assembly():
    records = [Record("a", "G" * 1000), Record("b", "A" * 500), Record("c", "C" * 200)]
    st = summarise(records)
    assert st.n_contigs == 3
    assert st.total_bp == 1700
    assert st.longest_bp == 1000
    assert st.shortest_bp == 200
    assert st.n50 == 1000
    assert st.l50 == 1


def test_summarise_copes_with_nothing_at_all():
    st = summarise([])
    assert st.n_contigs == 0 and st.total_bp == 0 and st.n50 == 0


def test_summary_lines_are_readable_sentences():
    lines = summarise([Record("a", "ACGT" * 100)]).summary_lines()
    assert any("N50" in line for line in lines)
    assert all(isinstance(line, str) and line for line in lines)


def test_read_length_summary_reports_read_n50():
    reads = ["A" * 10_000] + ["A" * 100] * 50
    info = read_length_summary(reads)
    assert info["n_reads"] == 51
    assert info["longest_bp"] == 10_000
    # The mean is dragged down by the crowd of short reads. N50 is not, which is
    # the reason it is the number worth looking at.
    assert info["n50_bp"] == 10_000
    assert info["mean_bp"] < 1_000


def test_read_length_summary_handles_no_reads():
    assert read_length_summary([])["n_reads"] == 0


def test_coverage_estimate_is_just_a_ratio():
    assert coverage_estimate(3_000_000, 100_000) == pytest.approx(30.0)
    assert coverage_estimate(100, 0) == 0.0
