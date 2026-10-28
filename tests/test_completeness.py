"""Tests for the optional BUSCO step.

BUSCO itself is not installed here and is not in environment.yml, so these tests
do not run it. What they do test is the part that can silently be wrong: reading
its summary. A parser that misses a field and returns a clean zero is worse than
one that crashes, because 0.0% complete looks like a real answer. So the fixtures
below are real BUSCO output in the shapes it actually writes, and the tests check
that every number comes back.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from nagasm import completeness

DATA = Path(__file__).parent / "data"


# A BUSCO 5 JSON summary, trimmed to the fields that matter.
JSON_SUMMARY = json.dumps(
    {
        "parameters": {"in": "nuclear.fasta", "mode": "genome"},
        "lineage_dataset": {"name": "eukaryota_odb10", "creation_date": "2024-01-08"},
        "results": {
            "one_line_summary": "C:93.3%[S:90.2%,D:3.1%],F:2.4%,M:4.3%,n:255",
            "Complete": 238,
            "Single copy": 230,
            "Multi copy": 8,
            "Fragmented": 6,
            "Missing": 11,
            "n_markers": 255,
            "Complete percentage": 93.3,
            "Single copy percentage": 90.2,
            "Multi copy percentage": 3.1,
            "Fragmented percentage": 2.4,
            "Missing percentage": 4.3,
        },
    }
)

# The plain text summary, in the layout BUSCO has used for years.
TEXT_SUMMARY = """# BUSCO version is: 5.8.2
# The lineage dataset is: bacteria_odb12 (Creation date: 2024-01-08, number of genomes: 4085, number of BUSCOs: 124)
# Summarized benchmarking in BUSCO notation for file nuclear.fasta
# BUSCO was run in mode: euk_genome_met

	***** Results: *****

	C:99.2%[S:98.4%,D:0.8%],F:0.0%,M:0.8%,n:124
	123	Complete BUSCOs (C)
	122	Complete and single-copy BUSCOs (S)
	1	Complete and duplicated BUSCOs (D)
	0	Fragmented BUSCOs (F)
	1	Missing BUSCOs (M)
	124	Total BUSCO groups searched
"""


def test_the_json_summary_is_read_in_full():
    got = completeness.parse_summary(JSON_SUMMARY)
    assert got["n_markers"] == 255
    assert got["complete_pct"] == pytest.approx(93.3)
    assert got["single_pct"] == pytest.approx(90.2)
    assert got["duplicated_pct"] == pytest.approx(3.1)
    assert got["fragmented_pct"] == pytest.approx(2.4)
    assert got["missing_pct"] == pytest.approx(4.3)
    assert got["complete_n"] == 238
    assert got["single_n"] == 230
    assert got["duplicated_n"] == 8
    assert got["fragmented_n"] == 6
    assert got["missing_n"] == 11


def test_the_text_summary_is_read_in_full():
    got = completeness.parse_summary(TEXT_SUMMARY)
    assert got["n_markers"] == 124
    assert got["complete_pct"] == pytest.approx(99.2)
    assert got["complete_n"] == 123
    assert got["single_n"] == 122
    assert got["duplicated_n"] == 1
    assert got["fragmented_n"] == 0
    assert got["missing_n"] == 1


def test_the_counts_are_recovered_when_only_the_one_line_is_there():
    # Some wrappers keep only that line. Guessing the counts back from the
    # shares is approximate, but a rounded count beats a zero that reads like a
    # measurement.
    got = completeness.parse_summary("C:80.0%[S:76.0%,D:4.0%],F:5.0%,M:15.0%,n:200")
    assert got["n_markers"] == 200
    assert got["complete_n"] == 160
    assert got["duplicated_n"] == 8
    assert got["missing_n"] == 30


def test_renamed_json_fields_still_parse():
    # BUSCO has moved these field names more than once, so the parser must fall
    # back to the one line summary rather than return zeros.
    payload = json.dumps(
        {"results": {"one_line_summary": "C:50.0%[S:49.0%,D:1.0%],F:10.0%,M:40.0%,n:100"}}
    )
    got = completeness.parse_summary(payload)
    assert got["complete_pct"] == pytest.approx(50.0)
    assert got["n_markers"] == 100


def test_a_summary_with_nothing_in_it_returns_nothing_not_zeros():
    # An empty dict means "could not read it". A dict full of zeros would be
    # reported as a real score of 0%, which is a lie.
    assert completeness.parse_summary("this is not a busco summary") == {}


def test_broken_json_falls_back_to_scanning_the_text():
    got = completeness.parse_summary('{"results": {"one_line_sum  C:70.0%[S:69.0%,D:1.0%],F:5.0%,M:25.0%,n:50')
    assert got["complete_pct"] == pytest.approx(70.0)
    assert got["n_markers"] == 50


# --------------------------------------------------------------- the summary file


def test_the_summary_file_is_found_wherever_busco_put_it(tmp_path):
    deep = tmp_path / "nuclear" / "run_eukaryota_odb10"
    deep.mkdir(parents=True)
    (deep / "short_summary.specific.eukaryota_odb10.nuclear.txt").write_text(TEXT_SUMMARY)
    found = completeness.find_summary(tmp_path)
    assert found is not None and found.name.startswith("short_summary")


def test_the_json_summary_is_preferred_over_the_text_one(tmp_path):
    folder = tmp_path / "nuclear"
    folder.mkdir()
    (folder / "short_summary.specific.eukaryota_odb10.nuclear.txt").write_text(TEXT_SUMMARY)
    (folder / "short_summary.specific.eukaryota_odb10.nuclear.json").write_text(JSON_SUMMARY)
    found = completeness.find_summary(tmp_path)
    assert found.suffix == ".json", "the JSON is the more reliable of the two"


def test_no_summary_anywhere_is_reported_as_none(tmp_path):
    assert completeness.find_summary(tmp_path) is None
    assert completeness.find_summary(tmp_path / "never_existed") is None


# ------------------------------------------------------------------ what it means


def test_a_small_lineage_is_called_coarse():
    result = completeness.Completeness(
        lineage="eukaryota_odb12", scored="nuclear.fasta", n_markers=255, complete_pct=95.0
    )
    notes = " ".join(completeness.interpret(result)).lower()
    assert "coarse" in notes, "a 255 marker set must not be presented as proof"


def test_scoring_the_nuclear_contigs_alone_is_stated():
    result = completeness.Completeness(
        lineage="eukaryota_odb12", scored="nuclear.fasta", n_markers=255, complete_pct=95.0
    )
    notes = " ".join(completeness.interpret(result)).lower()
    assert "nuclear contigs alone" in notes


def test_a_bacterial_run_is_not_told_about_its_nucleus():
    # The compartment is called nuclear inside the pipeline, because that is what
    # the depth triage separates. The report is read by people, and should not
    # hand a bacterium a nucleus.
    result = completeness.Completeness(
        lineage="bacteria_odb12", scored="nuclear.fasta", n_markers=124, complete_pct=99.2
    )
    notes = " ".join(completeness.interpret(result))
    assert "nuclear" not in notes.lower()
    assert "chromosomal" in notes


def test_a_close_small_lineage_is_not_dismissed_as_coarse():
    # 124 bacterial markers on a 5 Mb genome is not the same evidence as 255
    # markers shared across all eukaryotes, and calling both coarse would be lazy.
    close = completeness.Completeness(
        lineage="bacteria_odb12", scored="nuclear.fasta", n_markers=124, complete_pct=99.2
    )
    far = completeness.Completeness(
        lineage="eukaryota_odb12", scored="nuclear.fasta", n_markers=255, complete_pct=99.2
    )
    assert "coarse" not in " ".join(completeness.interpret(close)).lower()
    assert "coarse" in " ".join(completeness.interpret(far)).lower()


def test_heavy_duplication_is_explained_rather_than_condemned():
    result = completeness.Completeness(
        lineage="eukaryota_odb12",
        scored="nuclear.fasta",
        n_markers=255,
        complete_pct=95.0,
        duplicated_pct=22.0,
    )
    notes = " ".join(completeness.interpret(result))
    assert "22.0%" in notes
    assert "assembled twice" in notes


def test_a_low_score_says_to_check_the_lineage_first():
    result = completeness.Completeness(
        lineage="metazoa_odb12", scored="nuclear.fasta", n_markers=954, complete_pct=12.0
    )
    notes = " ".join(completeness.interpret(result)).lower()
    assert "poor" in notes and "match" in notes


def test_fragmentation_points_at_polishing():
    result = completeness.Completeness(
        lineage="eukaryota_odb12",
        scored="nuclear.fasta",
        n_markers=255,
        complete_pct=70.0,
        fragmented_pct=18.0,
    )
    notes = " ".join(completeness.interpret(result)).lower()
    assert "polish" in notes


def test_n50_is_always_mentioned_alongside():
    # The whole reason this module exists is that N50 and completeness can move
    # in opposite directions. Every report must say so.
    for pct in (5.0, 75.0, 99.0):
        result = completeness.Completeness(
            lineage="eukaryota_odb12", scored="nuclear.fasta", n_markers=255, complete_pct=pct
        )
        assert "N50" in " ".join(completeness.interpret(result))


def test_summary_lines_keep_the_counts_next_to_the_shares():
    result = completeness.Completeness(
        lineage="eukaryota_odb12",
        scored="nuclear.fasta",
        n_markers=255,
        complete_pct=93.3,
        complete_n=238,
    )
    text = "\n".join(result.summary_lines())
    assert "93.3%" in text and "238" in text and "255" in text


# ------------------------------------------------------------------- running it


def test_a_missing_input_file_is_caught_before_busco_is_looked_for(tmp_path):
    with pytest.raises(FileNotFoundError):
        completeness.run_busco(tmp_path / "nope.fasta", tmp_path / "out")


def test_no_busco_installed_says_how_to_install_it(tmp_path, monkeypatch):
    fasta = tmp_path / "nuclear.fasta"
    fasta.write_text(">c\nACGT\n")
    monkeypatch.setattr(completeness, "have_busco", lambda: False)
    with pytest.raises(completeness.BuscoMissing) as exc:
        completeness.run_busco(fasta, tmp_path / "out")
    assert "environment-busco.yml" in str(exc.value)


def test_the_suggested_lineages_all_explain_themselves():
    assert "eukaryota_odb12" in completeness.SUGGESTED_LINEAGES
    assert "bacteria_odb12" in completeness.SUGGESTED_LINEAGES
    for name, what in completeness.SUGGESTED_LINEAGES.items():
        assert what and len(what) > 5, f"{name} has no explanation"


def test_the_pipeline_skips_the_check_unless_a_lineage_is_asked_for(tmp_path):
    # The default must not reach for BUSCO at all, since it is not installed and
    # downloads a dataset when it is.
    from nagasm.pipeline import Settings, _score_completeness

    payload: dict = {}
    _score_completeness(Settings(), tmp_path, {}, payload, lambda msg: None)
    assert payload == {}, "the completeness check ran when nobody asked for it"


def test_asking_for_a_lineage_with_no_nuclear_contigs_is_a_note_not_a_crash(tmp_path):
    from nagasm.pipeline import Settings, _score_completeness

    payload: dict = {}
    messages: list[str] = []
    settings = Settings(busco="eukaryota_odb12")
    _score_completeness(settings, tmp_path, {}, payload, messages.append)
    assert "error" in payload["completeness"]
    assert messages, "the skip was silent"


def test_a_failing_busco_never_takes_the_run_down(tmp_path, monkeypatch):
    from nagasm.pipeline import Settings, _score_completeness

    nuclear = tmp_path / "nuclear.fasta"
    nuclear.write_text(">c\n" + "ACGT" * 100 + "\n")
    monkeypatch.setattr(
        completeness, "run_busco", lambda *a, **k: (_ for _ in ()).throw(RuntimeError("boom"))
    )
    payload: dict = {}
    settings = Settings(busco="eukaryota_odb12")
    _score_completeness(settings, tmp_path, {"nuclear": nuclear}, payload, lambda m: None)
    assert payload["completeness"]["error"] == "boom"


def test_the_report_states_which_file_was_scored(tmp_path):
    from nagasm import report

    payload = {
        "started": "today",
        "completeness": completeness.Completeness(
            lineage="eukaryota_odb12",
            scored="nuclear.fasta",
            n_markers=255,
            complete_pct=93.3,
            complete_n=238,
            one_line="C:93.3%[S:90.2%,D:3.1%],F:2.4%,M:4.3%,n:255",
            notes=["This score is for the nuclear contigs alone"],
        ).as_dict(),
    }
    text = report.write_report(payload, tmp_path / "r.md").read_text()
    assert "Gene content" in text
    assert "nuclear.fasta" in text
    assert "eukaryota_odb12" in text
    assert "93.3%" in text
    assert "nuclear contigs alone" in text


def test_a_failed_check_is_reported_as_failed_not_as_a_zero_score(tmp_path):
    from nagasm import report

    payload = {
        "started": "today",
        "completeness": {"error": "busco is not on PATH", "lineage": "eukaryota_odb12"},
    }
    text = report.write_report(payload, tmp_path / "r.md").read_text()
    assert "did not finish" in text
    assert "0.0%" not in text, "a failure must never be printed as a score"


# ------------------------------------------- against real BUSCO 6.0.0 output
#
# These two files came out of an actual run: BUSCO 6.0.0 scoring the E. coli
# assembly from the worked example against bacteria_odb12. They are here because
# everything above was written against output shapes taken from documentation,
# and documentation is not the same as the bytes a program writes.


def test_the_real_busco_6_json_parses_exactly():
    got = completeness.parse_summary((DATA / "busco6_short_summary.json").read_text())
    assert got["n_markers"] == 116
    assert got["complete_pct"] == pytest.approx(96.6)
    assert got["single_pct"] == pytest.approx(96.6)
    assert got["duplicated_pct"] == pytest.approx(0.0)
    assert got["fragmented_pct"] == pytest.approx(2.6)
    assert got["missing_pct"] == pytest.approx(0.9)
    assert got["complete_n"] == 112
    assert got["single_n"] == 112
    assert got["duplicated_n"] == 0
    assert got["fragmented_n"] == 3
    assert got["missing_n"] == 1


def test_the_real_busco_6_text_parses_exactly():
    got = completeness.parse_summary((DATA / "busco6_short_summary.txt").read_text())
    assert got["n_markers"] == 116
    assert got["complete_pct"] == pytest.approx(96.6)
    assert got["complete_n"] == 112
    assert got["single_n"] == 112
    assert got["duplicated_n"] == 0
    assert got["fragmented_n"] == 3
    assert got["missing_n"] == 1


def test_both_real_forms_agree_with_each_other():
    # The JSON and the text are two renderings of one run, so any field present
    # in both has to match. This is what catches a parser that reads the right
    # number from one and the wrong number from the other.
    a = completeness.parse_summary((DATA / "busco6_short_summary.json").read_text())
    b = completeness.parse_summary((DATA / "busco6_short_summary.txt").read_text())
    for key in set(a) & set(b) - {"one_line"}:
        assert a[key] == pytest.approx(b[key]), f"{key} differs: {a[key]} vs {b[key]}"


def test_a_real_zero_is_a_zero_and_not_a_missing_field():
    # This run had no duplicated markers. A parser that treats 0 as absent and
    # falls back to a guess would quietly invent a number here.
    got = completeness.parse_summary((DATA / "busco6_short_summary.json").read_text())
    assert got["duplicated_n"] == 0
    assert got["duplicated_pct"] == 0.0


# ------------------------------------------------------- the lineage name trap


def test_an_orthodb_release_suffix_is_stripped():
    # busco --list-datasets prints bacteria_odb12.2, but BUSCO builds the
    # download URL from bacteria_odb12 and then checks the hash listed against
    # bacteria_odb12.2. The file is fine, the comparison is not, and it fails
    # with "md5 hash is incorrect" after fetching 58 MB.
    name, note = completeness._normalise_lineage("bacteria_odb12.2")
    assert name == "bacteria_odb12"
    assert note and "md5" in note


def test_a_plain_lineage_name_is_left_alone():
    assert completeness._normalise_lineage("bacteria_odb12") == ("bacteria_odb12", None)
    assert completeness._normalise_lineage("eukaryota_odb10") == ("eukaryota_odb10", None)


def test_a_path_to_a_lineage_is_not_rewritten():
    name, note = completeness._normalise_lineage("/data/lineages/bacteria_odb12.2")
    assert name == "/data/lineages/bacteria_odb12.2"
    assert note is None


def test_the_suggested_lineages_are_all_names_busco_6_accepts():
    # Every one of these was checked against busco --list-datasets on 6.0.0.
    # None may carry a release suffix, since that is the failure above.
    for name in completeness.SUGGESTED_LINEAGES:
        assert completeness._normalise_lineage(name) == (name, None), name
    assert completeness.DEFAULT_LINEAGE in completeness.SUGGESTED_LINEAGES


# --------------------------------------------------- explaining real failures
#
# Each string below is from a log this module actually produced while it was
# being got working. BUSCO's own message points somewhere unhelpful in every
# case.


def test_the_md5_failure_is_explained_as_a_naming_problem():
    log = (
        "ERROR:\tmd5 hash is incorrect: a00b520c960409e7ffa01252e09cb733 while "
        "e21960d275e09687a86aace486f2886a expected\n"
        "INFO:\tdeleting corrupted file busco_downloads/lineages/bacteria_odb12.2.tar.gz\n"
    )
    hint = completeness.diagnose_failure(log)
    assert hint and "suffix" in hint
    assert "corrupt" not in hint.lower().replace("corrupted", "")


def test_the_missing_java_failure_names_java():
    log = "The operation couldn't be completed. Unable to locate a Java Runtime."
    hint = completeness.diagnose_failure(log)
    assert hint and "Java" in hint and "BBMap" in hint


def test_the_wrong_python_failure_names_the_path_order():
    hint = completeness.diagnose_failure("No module named 'busco'")
    assert hint and "PATH" in hint


def test_an_unrecognised_failure_gets_no_invented_explanation():
    assert completeness.diagnose_failure("something nobody has seen before") is None
    assert completeness.diagnose_failure("") is None
