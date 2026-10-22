"""Cleaning up a draft with a real aligner.

The assembler in `circular` lays reads out using anchor chains, which is enough
to get the structure and the length right but leaves the draft carrying whatever
error rate the reads had. Turning that into an accurate sequence needs real base
by base alignment, and that is a solved problem, so we hand it to minimap2 and
racon rather than writing a third rate version of it here.

This is the same split miniasm and racon use, and it is a deliberate one. The
parts of this package that are worth reading are the ones that do something a
general tool does not. Consensus polishing is not one of those parts.

## Polishing a circle without leaving a scar

A circle written out as a line has a join in it, and that join is the one place
no read can map across, because a read that spans it looks like two pieces at
opposite ends of the draft. Polish it as it stands and every position improves
except the few hundred bases either side of the join, which stay as rough as the
reads were.

The fix is to polish a doubled copy of the sequence. Every position of the circle
then sits somewhere in the middle of the doubled draft with reads spanning it
properly, and the real sequence is taken from the middle afterwards. The join
moves to a place where it does no harm instead of being left in the answer.
"""

from __future__ import annotations

import shutil
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path

from .seqio import Record, read_fasta, write_fasta, write_fastq


class ToolMissing(RuntimeError):
    """A command line tool this step needs is not on PATH."""


def have(tool: str) -> bool:
    return shutil.which(tool) is not None


def require(*tools: str) -> None:
    missing = [t for t in tools if not have(t)]
    if missing:
        raise ToolMissing(
            f"these tools are needed but not installed: {', '.join(missing)}. "
            f"They come with the environment in environment.yml."
        )


@dataclass
class PolishResult:
    sequence: str
    rounds: int
    changed_bases: int
    notes: list[str]


def _run(cmd: list[str], stdout_path: Path | None = None, log: Path | None = None) -> None:
    """Run a command, sending output where asked and raising with its own words."""
    out = open(stdout_path, "wb") if stdout_path else subprocess.DEVNULL
    err = open(log, "wb") if log else subprocess.PIPE
    try:
        proc = subprocess.run(cmd, stdout=out, stderr=err)
        if proc.returncode != 0:
            detail = ""
            if log and Path(log).exists():
                detail = Path(log).read_text()[-2000:]
            elif proc.stderr:
                detail = proc.stderr.decode()[-2000:]
            raise RuntimeError(f"{cmd[0]} failed:\n{detail}")
    finally:
        if stdout_path:
            out.close()
        if log:
            err.close()


def _preset_for(read_type: str) -> str:
    """The minimap2 preset that matches the reads."""
    return {
        "ont": "map-ont",
        "nanopore": "map-ont",
        "hifi": "map-hifi",
        "pacbio": "map-pb",
        "illumina": "sr",
        "short": "sr",
    }.get(read_type.lower(), "map-ont")


def polish(
    draft: str,
    reads: list[str],
    read_type: str = "ont",
    rounds: int = 2,
    threads: int = 4,
    circular: bool = False,
    workdir: Path | None = None,
) -> PolishResult:
    """Polish a draft with racon, mapping the reads back with minimap2.

    `circular` doubles the draft before polishing and takes the middle back out
    afterwards, so the join gets the same treatment as everywhere else.

    Two rounds is the usual advice. The second picks up places the first one's
    changes made mappable, and a third rarely earns its time.
    """
    require("minimap2", "racon")
    notes: list[str] = []
    if not draft or not reads:
        return PolishResult(draft, 0, 0, ["nothing to polish"])

    original_length = len(draft)
    working = draft + draft if circular else draft
    if circular:
        notes.append(
            "polished a doubled copy so that reads could span the join, then took "
            "the middle back out"
        )

    tmp = Path(workdir) if workdir else Path(tempfile.mkdtemp(prefix="nagasm_polish_"))
    tmp.mkdir(parents=True, exist_ok=True)

    reads_path = tmp / "reads.fastq"
    write_fastq(reads_path, [Record(f"read{i}", r) for i, r in enumerate(reads)])

    preset = _preset_for(read_type)
    current = working
    done = 0
    for round_number in range(1, rounds + 1):
        draft_path = tmp / f"draft{round_number}.fasta"
        write_fasta(draft_path, [Record("draft", current)], wrap=0)
        paf = tmp / f"map{round_number}.paf"
        _run(
            ["minimap2", "-x", preset, "-t", str(threads), str(draft_path), str(reads_path)],
            stdout_path=paf,
            log=tmp / f"minimap2_{round_number}.log",
        )
        if paf.stat().st_size == 0:
            notes.append(
                f"round {round_number}: no reads mapped back to the draft, so polishing "
                f"stopped here"
            )
            break
        out = tmp / f"polished{round_number}.fasta"
        try:
            _run(
                [
                    "racon",
                    "-t",
                    str(threads),
                    str(reads_path),
                    str(paf),
                    str(draft_path),
                ],
                stdout_path=out,
                log=tmp / f"racon_{round_number}.log",
            )
        except RuntimeError as exc:
            notes.append(f"round {round_number}: racon stopped early ({exc})")
            break
        records = list(read_fasta(out))
        if not records:
            notes.append(f"round {round_number}: racon returned nothing, keeping the draft")
            break
        # racon can drop a piece it could not support. Keep the longest.
        current = max((r.seq for r in records), key=len)
        done = round_number

    if circular:
        # Take a window of the original length from the middle, so the join that
        # was at the ends ends up safely inside.
        if len(current) >= original_length:
            start = (len(current) - original_length) // 2
            current = current[start : start + original_length]
        else:
            notes.append(
                "the polished doubled draft came back shorter than expected, so it "
                "was kept whole rather than cut"
            )

    changed = sum(1 for a, b in zip(draft, current) if a != b) + abs(len(draft) - len(current))
    return PolishResult(current, done, changed, notes)
