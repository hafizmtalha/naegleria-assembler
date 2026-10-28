"""Gene content, checked with BUSCO. Optional, and off unless asked for.

Everything else in this package measures the shape of an assembly: how long the
pieces are, how deep the reads sit on them, whether a circle closes, whether the
reads that support an element really lie inside it. None of that looks at genes.
An assembly can be long, well supported and badly incomplete at the same time,
and the usual way to notice is to look for genes that every member of a group is
expected to carry. That is what BUSCO does, so this module hands the job to it
rather than reinventing it.

Three things about running it here are worth saying plainly, because getting any
of them wrong turns the score into a number that means nothing.

**It scores the nuclear compartment, not the whole output.** The high copy
elements are not supposed to carry single copy nuclear genes. If you point BUSCO
at the elements and the nuclear contigs together, any marker that happens to sit
on both comes back as duplicated, and the duplication is not a duplication, it is
two different kinds of DNA in one file. The pipeline only ever scores
`nuclear.fasta`, and the report says so next to the score.

**The lineage is a compromise on this organism.** Naegleria sits in
Heterolobosea, inside Discoba, and BUSCO has no set anywhere near that close.
`eukaryota_odb10` is the honest choice, about 255 markers shared across all
eukaryotes. It is coarse. A good score on 255 broad markers is weak evidence that
a 30 Mb genome is complete, and the report phrases it that way rather than
calling the assembly finished. For the bacterial test case `bacteria_odb10` is
the right set and is much more informative, since it is both larger and far
closer to the organism.

**It needs to download its lineage the first time.** That is the one thing in
this pipeline that reaches out to the network, which is why it is off by default.
Pass `--busco-download-path` at a folder that already holds the lineage and it
runs offline.

The failure mode this is actually good at catching is the one N50 hides. An
assembler that collapses a repeat family produces fewer, longer contigs and a
better N50 while quietly losing sequence. Contiguity goes up, completeness goes
down. Reading the two together is the point.
"""

from __future__ import annotations

import json
import re
import shutil
import subprocess
from dataclasses import asdict, dataclass, field
from pathlib import Path

# The lineages worth suggesting, with what they are for. Not a complete list,
# BUSCO has hundreds, and `busco --list-datasets` is the only authority on what
# your install can actually fetch.
#
# The suffix is not cosmetic, and it cost an hour. Two traps sit here, both found
# by installing BUSCO 6.0.0 and watching it fail.
#
# First, BUSCO 6 dropped every odb10 set for cellular life and ships only odb12.
# Only the viral sets kept the odb10 name. So `eukaryota_odb10`, which is what
# most papers quote and what this module first defaulted to, is not a name BUSCO
# 6 will accept.
#
# Second, and worse, `busco --list-datasets` prints names ending `_odb12.2`, but
# those are the OrthoDB release label rather than the dataset name. Ask for
# `bacteria_odb12.2` and BUSCO builds the download URL for `bacteria_odb12`,
# fetches the right 58 MB file, then checks it against the md5 on the
# `bacteria_odb12.2` row of its own manifest, which is a different hash. It fails
# with "md5 hash is incorrect" and deletes the file, which reads like a corrupt
# download and is not one. The name to pass is `bacteria_odb12`, with no `.2`.
# See _normalise_lineage below, which strips it for you.
#
# The marker counts are the ones BUSCO itself reports, not estimates. On BUSCO 5
# you want the odb10 spellings and different counts: eukaryota_odb10 held 255
# markers where eukaryota_odb12 holds 125. The two generations score the same
# assembly differently, so whichever you use, say which one in the methods.
DEFAULT_LINEAGE = "eukaryota_odb12"

SUGGESTED_LINEAGES = {
    "eukaryota_odb12": "all eukaryotes, 125 markers. The honest choice for "
    "Naegleria, and a very coarse one. BUSCO 5 called this eukaryota_odb10 and "
    "had 255 markers in it",
    "bacteria_odb12": "all bacteria, 116 markers. The right set for the "
    "bacterial test case, and a close match to the organism",
    "archaea_odb12": "all archaea, 195 markers",
    "euglenozoa_odb12": "Euglenozoa, 2,579 markers. The other main branch of "
    "Discoba, so it is Naegleria's nearest relative in the list and still far "
    "too specific to use on it. Here to be explicit that nothing closer exists",
    "metazoa_odb12": "animals, 932 markers. Wrong for Naegleria, and a useful "
    "illustration of what a badly matched lineage scores",
    "fungi_odb12": "fungi, 1,019 markers. A common point of comparison, not a "
    "fit for this organism",
}

# The OrthoDB release label that `busco --list-datasets` prints but BUSCO will
# not accept as a dataset name. See the note above.
_RELEASE_LABEL = re.compile(r"^(?P<name>.+_odb\d+)\.\d+$")


def _normalise_lineage(lineage: str) -> tuple[str, str | None]:
    """Turn a name copied out of `--list-datasets` into one BUSCO will accept.

    Returns the name to use and, when something was changed, a note explaining
    what. Copying `bacteria_odb12.2` straight out of the dataset listing is the
    obvious thing to do and it fails in a way that looks like a corrupt download,
    so this fixes it quietly and says so rather than letting the run die twenty
    seconds and 58 MB later.

    A path is left alone, since an explicit directory is not a dataset name.
    """
    if "/" in lineage or Path(lineage).exists():
        return lineage, None
    found = _RELEASE_LABEL.match(lineage)
    if not found:
        return lineage, None
    name = found.group("name")
    return name, (
        f"asked for {lineage}, using {name}. The trailing number is the OrthoDB "
        "release that busco --list-datasets prints, not part of the dataset name, "
        "and passing it makes BUSCO fail its own md5 check"
    )

ONE_LINE = re.compile(
    r"C:\s*([\d.]+)%\s*\[\s*S:\s*([\d.]+)%\s*,\s*D:\s*([\d.]+)%\s*\]\s*,\s*"
    r"F:\s*([\d.]+)%\s*,\s*M:\s*([\d.]+)%\s*,\s*n:\s*(\d+)"
)


class BuscoMissing(RuntimeError):
    """BUSCO is not on PATH."""


@dataclass
class Completeness:
    """What BUSCO found, with the counts and the shares kept apart.

    The shares are what BUSCO prints and what people quote. The counts are what
    you need to say anything careful, because 95% of 255 markers and 95% of 2,000
    markers are not the same evidence.
    """

    lineage: str
    scored: str  # which file was scored, so a score is never read as whole-assembly
    n_markers: int = 0
    complete_pct: float = 0.0
    single_pct: float = 0.0
    duplicated_pct: float = 0.0
    fragmented_pct: float = 0.0
    missing_pct: float = 0.0
    complete_n: int = 0
    single_n: int = 0
    duplicated_n: int = 0
    fragmented_n: int = 0
    missing_n: int = 0
    one_line: str = ""
    outdir: str = ""
    notes: list[str] = field(default_factory=list)

    def as_dict(self) -> dict:
        return asdict(self)

    def summary_lines(self) -> list[str]:
        return [
            f"Lineage: {self.lineage} ({self.n_markers:,} markers)",
            f"Scored: {self.scored}",
            f"Complete: {self.complete_pct:.1f}% ({self.complete_n:,})",
            f"  of those, single copy: {self.single_pct:.1f}% ({self.single_n:,})",
            f"  of those, duplicated: {self.duplicated_pct:.1f}% ({self.duplicated_n:,})",
            f"Fragmented: {self.fragmented_pct:.1f}% ({self.fragmented_n:,})",
            f"Missing: {self.missing_pct:.1f}% ({self.missing_n:,})",
        ]


def have_busco() -> bool:
    return shutil.which("busco") is not None


def busco_version() -> str | None:
    """The installed version, or None if BUSCO is not there or will not answer."""
    if not have_busco():
        return None
    try:
        out = subprocess.run(
            ["busco", "--version"], capture_output=True, text=True, timeout=60
        )
    except (OSError, subprocess.SubprocessError):
        return None
    # Some versions print the version to stdout, some to stderr.
    text = (out.stdout or "") + (out.stderr or "")
    found = re.search(r"(\d+\.\d+(?:\.\d+)?)", text)
    return found.group(1) if found else None


def _as_int(value) -> int:
    """A count out of whatever BUSCO put in the field, including None or a string."""
    try:
        return int(float(value))
    except (TypeError, ValueError):
        return 0


def _as_float(value) -> float:
    try:
        return float(str(value).rstrip("%"))
    except (TypeError, ValueError):
        return 0.0


def _first(mapping: dict, *keys):
    """The first key present, since the field names move between BUSCO versions."""
    for key in keys:
        if key in mapping and mapping[key] is not None:
            return mapping[key]
    return None


def parse_summary(text: str) -> dict:
    """Pull the numbers out of a BUSCO short summary, JSON or plain text.

    Both forms carry the same one line summary, in the shape
    ``C:95.0%[S:94.0%,D:1.0%],F:2.0%,M:3.0%,n:255``. That line has not changed
    across BUSCO versions, so it is the primary source here. The JSON field names
    have changed, more than once, so they are read with fallbacks and only used
    to fill in the counts the one line summary does not carry.
    """
    numbers: dict = {}

    payload = None
    stripped = text.lstrip()
    if stripped.startswith("{"):
        try:
            payload = json.loads(text)
        except json.JSONDecodeError:
            payload = None

    results = {}
    if isinstance(payload, dict):
        results = payload.get("results") or payload

    one_line = ""
    if isinstance(results, dict):
        one_line = str(_first(results, "one_line_summary", "one line summary") or "")
    if not one_line:
        for line in text.splitlines():
            if ONE_LINE.search(line):
                one_line = line.strip()
                break

    found = ONE_LINE.search(one_line) or ONE_LINE.search(text)
    if found:
        numbers.update(
            complete_pct=float(found.group(1)),
            single_pct=float(found.group(2)),
            duplicated_pct=float(found.group(3)),
            fragmented_pct=float(found.group(4)),
            missing_pct=float(found.group(5)),
            n_markers=int(found.group(6)),
            one_line=found.group(0),
        )

    if isinstance(results, dict) and results:
        counts = {
            "complete_n": ("Complete BUSCOs", "Complete"),
            "single_n": ("Single copy BUSCOs", "Single copy"),
            "duplicated_n": ("Multi copy BUSCOs", "Multi copy", "Duplicated"),
            "fragmented_n": ("Fragmented BUSCOs", "Fragmented"),
            "missing_n": ("Missing BUSCOs", "Missing"),
        }
        for field_name, keys in counts.items():
            value = _first(results, *keys)
            if value is not None:
                numbers[field_name] = _as_int(value)
        markers = _first(results, "n_markers", "n")
        if markers is not None and not numbers.get("n_markers"):
            numbers["n_markers"] = _as_int(markers)
        for field_name, key in (
            ("complete_pct", "Complete percentage"),
            ("single_pct", "Single copy percentage"),
            ("duplicated_pct", "Multi copy percentage"),
            ("fragmented_pct", "Fragmented percentage"),
            ("missing_pct", "Missing percentage"),
        ):
            if key in results and field_name not in numbers:
                numbers[field_name] = _as_float(results[key])

    # The plain text summary carries the counts as lines like
    # "        242     Complete BUSCOs (C)". Read them if the JSON was absent.
    if not numbers.get("complete_n"):
        patterns = {
            "complete_n": r"(\d+)\s+Complete BUSCOs",
            "single_n": r"(\d+)\s+Complete and single-copy BUSCOs",
            "duplicated_n": r"(\d+)\s+Complete and duplicated BUSCOs",
            "fragmented_n": r"(\d+)\s+Fragmented BUSCOs",
            "missing_n": r"(\d+)\s+Missing BUSCOs",
            "n_markers": r"(\d+)\s+Total BUSCO groups searched",
        }
        for field_name, pattern in patterns.items():
            hit = re.search(pattern, text)
            if hit:
                numbers.setdefault(field_name, int(hit.group(1)))

    # Last resort for the counts: the shares and the marker total give them back,
    # and a rounded count is better than a zero that reads like a real answer.
    total = numbers.get("n_markers", 0)
    if total and not numbers.get("complete_n"):
        for count_field, pct_field in (
            ("complete_n", "complete_pct"),
            ("single_n", "single_pct"),
            ("duplicated_n", "duplicated_pct"),
            ("fragmented_n", "fragmented_pct"),
            ("missing_n", "missing_pct"),
        ):
            numbers[count_field] = round(total * numbers.get(pct_field, 0.0) / 100)

    return numbers


def find_summary(outdir: Path) -> Path | None:
    """The short summary BUSCO wrote, preferring the JSON.

    BUSCO puts it at ``<out_path>/<run name>/short_summary.specific.<lineage>.<run
    name>.json``, and has moved the exact path around between versions, so this
    looks rather than assumes.
    """
    outdir = Path(outdir)
    if not outdir.exists():
        return None
    for pattern in ("short_summary*.json", "short_summary*.txt", "short_summary*"):
        found = sorted(outdir.rglob(pattern))
        if found:
            return found[0]
    return None


def diagnose_failure(log_text: str) -> str | None:
    """Turn a BUSCO failure log into the thing you actually need to change.

    Every entry here is a failure that happened while this module was being
    written, and each one has a log message that points somewhere unhelpful. The
    md5 one is the worst: it reads like a corrupt download, deletes the file, and
    is really a wrong name. Handing someone BUSCO's own last twenty lines and
    letting them work that out again is not good enough.
    """
    if not log_text:
        return None

    if "md5 hash is incorrect" in log_text:
        return (
            "BUSCO downloaded the lineage and then failed its own md5 check. This "
            "almost always means the lineage name carried an OrthoDB release "
            "suffix, like bacteria_odb12.2 instead of bacteria_odb12. BUSCO builds "
            "the URL from the name without the suffix but checks the hash listed "
            "against the name with it, so the file is fine and the comparison is "
            "not. Drop the trailing number"
        )

    if "Unable to locate a Java Runtime" in log_text or "java: command not found" in log_text:
        return (
            "BUSCO needs BBMap for its assembly statistics, and BBMap is Java. "
            "There is no Java on PATH. Install openjdk into the BUSCO "
            "environment, and note that some conda builds put java in "
            "lib/jvm/bin rather than bin, so activate the environment properly "
            "rather than only prepending its bin folder"
        )

    if "No module named 'busco'" in log_text:
        return (
            "BUSCO's launcher starts with #!/usr/bin/env python3, so it runs "
            "whichever python3 is first on PATH. Another environment is ahead of "
            "it and does not have BUSCO installed. Put the BUSCO environment's "
            "bin first. The nagasm command has an absolute interpreter path, so "
            "it keeps working either way"
        )

    if "unable to download or update all necessary files" in log_text.lower():
        return (
            "BUSCO could not fetch the lineage. If the machine is offline, "
            "download it once elsewhere and point --busco-download-path at the "
            "folder holding it"
        )

    if "not a valid option for" in log_text or "Unknown dataset" in log_text:
        return (
            "BUSCO does not recognise that lineage. Run busco --list-datasets to "
            "see what your version offers, and remember that BUSCO 6 dropped the "
            "odb10 sets for everything except viruses"
        )

    return None


def interpret(result: Completeness, scored_nuclear_only: bool = True) -> list[str]:
    """Plain sentences about what the score does and does not show.

    A bare percentage invites two opposite mistakes: reading a high score as a
    finished genome, and reading duplication as a problem when it is often the
    assembly telling you something true about the organism. These notes are the
    things worth saying before anyone quotes the number.
    """
    notes: list[str] = []

    # On a prokaryote there is no nucleus and no nuclear gene, so the wording has
    # to follow the lineage. The compartment is still called nuclear inside this
    # pipeline, because that is what the depth triage separates, but the report is
    # read by people and should not claim a nucleus the organism does not have.
    prokaryote = result.lineage.startswith(("bacteria", "archaea"))
    compartment = "chromosomal" if prokaryote else "nuclear"

    if scored_nuclear_only:
        notes.append(
            f"This score is for the {compartment} contigs alone. The high copy "
            "elements were left out on purpose, since they are not expected to "
            f"carry single copy {compartment} genes and including them would turn "
            "ordinary markers into false duplicates"
        )

    if result.n_markers and result.n_markers < 400:
        if prokaryote:
            # A small set that sits right on top of the organism is a different
            # thing from a small set stretched across all eukaryotes, and saying
            # both are "coarse" would be lazy.
            notes.append(
                f"The lineage has {result.n_markers:,} markers, which is a small set, "
                "but it is matched closely to the organism and spread over a genome "
                "of a few megabases, so it carries real weight here"
            )
        else:
            notes.append(
                f"The lineage has only {result.n_markers:,} markers, shared across all "
                "of a very wide group, so this is a coarse check. Spread over a genome "
                "this size, a high score is weak evidence of completeness rather than "
                "proof of it, while a low score is still worth taking seriously"
            )

    if result.complete_pct >= 90:
        notes.append(
            f"{result.complete_pct:.1f}% of the markers are present and complete, which "
            "says the assembly is not missing large tracts of conserved gene space"
        )
    elif result.complete_pct >= 70:
        notes.append(
            f"{result.complete_pct:.1f}% complete leaves a real gap. Check whether the "
            "missing markers are absent from the reads as well, or only from the "
            "assembly, since those have different causes"
        )
    elif result.complete_pct > 0:
        notes.append(
            f"{result.complete_pct:.1f}% complete is low. Either the lineage is a poor "
            "match for this organism, or a large share of the genome is not in the "
            "assembly. Worth settling which before reading anything else here"
        )

    if result.duplicated_pct >= 10:
        notes.append(
            f"{result.duplicated_pct:.1f}% of the markers turn up more than once. On a "
            "genome that has not been collapsed to one haplotype this is expected. If "
            "it is not that, it points at the same sequence assembled twice"
        )

    if result.fragmented_pct >= 10:
        notes.append(
            f"{result.fragmented_pct:.1f}% of the markers are fragmented, which usually "
            "means base level errors breaking the gene models rather than missing "
            "sequence. More polishing rounds are the first thing to try"
        )

    notes.append(
        "Read this next to N50, not instead of it. An assembler that collapses a "
        "repeat family gets a better N50 and a worse completeness at the same time, "
        "so the two numbers moving in opposite directions is itself the finding"
    )
    return notes


def run_busco(
    fasta: Path | str,
    outdir: Path | str,
    lineage: str = DEFAULT_LINEAGE,
    threads: int = 4,
    mode: str = "genome",
    download_path: Path | str | None = None,
    name: str = "nuclear",
    extra_args: list[str] | None = None,
    timeout: int | None = None,
) -> Completeness:
    """Run BUSCO on one FASTA and read the summary back.

    ``download_path`` points at a folder that already holds the lineage, and
    turns on BUSCO's offline mode. Without it BUSCO fetches the lineage the first
    time, which is the only step in this pipeline that touches the network.

    Raises BuscoMissing if BUSCO is not installed and RuntimeError if it ran and
    failed. Both are caught by the pipeline, which never lets this step take a
    finished assembly down with it.
    """
    fasta = Path(fasta)
    outdir = Path(outdir)
    if not fasta.exists():
        raise FileNotFoundError(f"nothing to score at {fasta}")
    if not have_busco():
        raise BuscoMissing(
            "busco is not on PATH. It is deliberately not in environment.yml, "
            "because it pulls in a large set of gene prediction tools that can "
            "fight with the rest of this environment. Install it on its own with "
            "environment-busco.yml, or with: mamba create -n busco -c conda-forge "
            "-c bioconda busco=6.0.0 sepp=4.5.5. The sepp pin matters, BUSCO's own "
            "guide warns that 4.5.6 picks the wrong lineage without saying so"
        )

    lineage, renamed = _normalise_lineage(str(lineage))
    extra_notes = [renamed] if renamed else []

    outdir.mkdir(parents=True, exist_ok=True)
    cmd = [
        "busco",
        "-i",
        str(fasta),
        "-o",
        name,
        "--out_path",
        str(outdir),
        "-m",
        mode,
        "-l",
        str(lineage),
        "-c",
        str(threads),
        "-f",  # overwrite, so a second run of the pipeline does not stall
    ]
    if download_path:
        cmd += ["--offline", "--download_path", str(download_path)]
    if extra_args:
        cmd += list(extra_args)

    log_path = outdir / "busco.log"
    with open(log_path, "wb") as handle:
        proc = subprocess.run(cmd, stdout=handle, stderr=subprocess.STDOUT, timeout=timeout)

    summary = find_summary(outdir / name) or find_summary(outdir)
    if proc.returncode != 0 and summary is None:
        log_text = log_path.read_text(errors="replace") if log_path.exists() else ""
        hint = diagnose_failure(log_text)
        message = f"busco failed, its log ends:\n{log_text[-2000:]}"
        if hint:
            message = f"busco failed: {hint}\n\nIts log ends:\n{log_text[-2000:]}"
        raise RuntimeError(message)
    if summary is None:
        raise RuntimeError(
            f"busco finished but wrote no short summary under {outdir}. Its log is "
            f"at {log_path}"
        )

    numbers = parse_summary(summary.read_text(errors="replace"))
    result = Completeness(
        lineage=str(lineage),
        scored=fasta.name,
        outdir=str(outdir / name),
        **{k: v for k, v in numbers.items() if k in Completeness.__dataclass_fields__},
    )
    result.notes = extra_notes + interpret(
        result, scored_nuclear_only=fasta.name == "nuclear.fasta"
    )
    return result
