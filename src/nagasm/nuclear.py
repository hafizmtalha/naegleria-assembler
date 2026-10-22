"""Assembling the nuclear compartment.

By the time this runs, triage has taken the high copy circles out, so what is
left is an ordinary genome at a roughly even depth. That is what mature
assemblers are built for, and in stage one this module hands the work to one of
them.

That is a deliberate choice rather than a gap. Flye and SPAdes represent many
years of work on repeat resolution and error handling, and a hand written
replacement would be slower and worse. The part of this package worth reading is
what happens before and after: deciding which reads belong here in the first
place, and giving each compartment the treatment it needs.

Stage two adds an assembler of our own, in `olc.py`, for people who want to read
one rather than use one. You can pick it with `engine="nagasm"`.

## Which tool for which reads

Three kinds of input are supported, and the right tool differs for each.

* Long reads on their own go to Flye, which was built for them.
* Short reads on their own go to SPAdes.
* Both together go to SPAdes in hybrid mode, where the long reads bridge repeats
  and the short reads supply the accuracy.

There is a fourth route worth knowing about, which is to assemble the long reads
with Flye and then polish with the short ones. On an AT rich genome read with
nanopore that often ends up better than hybrid SPAdes, because the errors left in
a long read assembly are mostly in homopolymers and short reads fix those
cleanly. It is available as `mode="long_then_polish"`.
"""

from __future__ import annotations

import shutil
import subprocess
from dataclasses import dataclass, field
from pathlib import Path

from .polish import ToolMissing, have, require
from .seqio import Record, read_fasta, write_fasta, write_fastq
from .stats import AssemblyStats, summarise


@dataclass
class NuclearAssembly:
    contigs: list[Record]
    engine: str
    mode: str
    stats: AssemblyStats | None = None
    notes: list[str] = field(default_factory=list)

    def __len__(self) -> int:
        return len(self.contigs)


def _run(cmd: list[str], log: Path) -> None:
    log.parent.mkdir(parents=True, exist_ok=True)
    with open(log, "wb") as handle:
        proc = subprocess.run(cmd, stdout=handle, stderr=subprocess.STDOUT)
    if proc.returncode != 0:
        tail = log.read_text(errors="replace")[-3000:]
        raise RuntimeError(f"{cmd[0]} failed. The end of its log says:\n{tail}")


def choose_mode(long_reads: list[str], short_reads: list[str]) -> str:
    """Work out which of the three input shapes we have."""
    if long_reads and short_reads:
        return "hybrid"
    if long_reads:
        return "long"
    if short_reads:
        return "short"
    raise ValueError("no reads given, so there is nothing to assemble")


def _flye_read_flag(read_type: str) -> str:
    """The Flye flag matching the platform and chemistry."""
    return {
        "ont": "--nano-raw",
        "nanopore": "--nano-raw",
        "ont-hq": "--nano-hq",
        "hifi": "--pacbio-hifi",
        "pacbio": "--pacbio-raw",
    }.get(read_type.lower(), "--nano-raw")


def assemble_with_flye(
    reads: list[str],
    outdir: Path,
    genome_size: int,
    read_type: str = "ont",
    threads: int = 4,
    keep_haplotypes: bool = False,
) -> list[Record]:
    """Run Flye on long reads.

    `genome_size` is a hint, and the estimate triage made for the nuclear
    compartment is a good one to pass. It does not have to be exact.

    `keep_haplotypes` is off by default. Naegleria has been reported to vary in
    ploidy, and on a genome like that keeping haplotypes apart can double parts of
    the assembly and make the total length look wrong. Turning it on is the right
    move if you are specifically looking for that variation.
    """
    require("flye")
    outdir.mkdir(parents=True, exist_ok=True)
    reads_path = outdir / "long_reads.fasta"
    write_fasta(reads_path, [Record(f"r{i}", r) for i, r in enumerate(reads)], wrap=0)

    cmd = [
        "flye",
        _flye_read_flag(read_type),
        str(reads_path),
        "--out-dir",
        str(outdir / "flye"),
        "--threads",
        str(threads),
    ]
    if genome_size:
        cmd += ["--genome-size", str(genome_size)]
    if keep_haplotypes:
        cmd += ["--keep-haplotypes"]

    _run(cmd, outdir / "flye.log")
    assembly = outdir / "flye" / "assembly.fasta"
    if not assembly.exists():
        raise RuntimeError(f"Flye finished but wrote no assembly at {assembly}")
    return list(read_fasta(assembly))


def assemble_with_spades(
    short_reads: list[str],
    outdir: Path,
    long_reads: list[str] | None = None,
    long_read_type: str = "ont",
    threads: int = 4,
    memory_gb: int = 12,
) -> list[Record]:
    """Run SPAdes on short reads, optionally with long reads to bridge repeats.

    The short reads are written as one interleaved file. If your pairs are in two
    files, interleave them before calling this, or pass them through the CLI,
    which handles it.
    """
    require("spades.py")
    outdir.mkdir(parents=True, exist_ok=True)
    short_path = outdir / "short_reads.fastq"
    write_fastq(short_path, [Record(f"r{i}", r) for i, r in enumerate(short_reads)])

    cmd = [
        "spades.py",
        "-s",
        str(short_path),
        "-o",
        str(outdir / "spades"),
        "--threads",
        str(threads),
        "--memory",
        str(memory_gb),
    ]
    if long_reads:
        long_path = outdir / "long_reads.fasta"
        write_fasta(
            long_path, [Record(f"l{i}", r) for i, r in enumerate(long_reads)], wrap=0
        )
        flag = "--pacbio" if "pac" in long_read_type.lower() or "hifi" in long_read_type.lower() else "--nanopore"
        cmd += [flag, str(long_path)]

    _run(cmd, outdir / "spades.log")
    for name in ("scaffolds.fasta", "contigs.fasta"):
        path = outdir / "spades" / name
        if path.exists():
            return list(read_fasta(path))
    raise RuntimeError("SPAdes finished but wrote neither scaffolds nor contigs")


def assemble_nuclear(
    long_reads: list[str] | None = None,
    short_reads: list[str] | None = None,
    outdir: Path = Path("results/nuclear"),
    genome_size: int = 0,
    read_type: str = "ont",
    mode: str | None = None,
    engine: str = "auto",
    threads: int = 4,
    memory_gb: int = 12,
    min_contig_bp: int = 500,
) -> NuclearAssembly:
    """Assemble the nuclear compartment with whichever route fits the reads.

    `engine` is "auto" to pick a proven assembler, or "nagasm" to use the one in
    this package. The second is there to be read and understood rather than to
    compete, and it says so in the notes.
    """
    long_reads = list(long_reads or [])
    short_reads = list(short_reads or [])
    mode = mode or choose_mode(long_reads, short_reads)
    outdir = Path(outdir)
    notes: list[str] = []

    if engine == "nagasm":
        from .olc import assemble_olc

        reads = long_reads or short_reads
        notes.append(
            "assembled with the assembler in this package. It is here to be read "
            "rather than to beat a mature tool, so expect a more broken up result "
            "than Flye or SPAdes would give"
        )
        contigs = assemble_olc(reads, min_contig_bp=min_contig_bp)
        used_engine = "nagasm"
    elif mode == "long":
        contigs = assemble_with_flye(
            long_reads, outdir, genome_size, read_type, threads
        )
        used_engine = "flye"
    elif mode == "short":
        contigs = assemble_with_spades(
            short_reads, outdir, threads=threads, memory_gb=memory_gb
        )
        used_engine = "spades"
    elif mode == "hybrid":
        contigs = assemble_with_spades(
            short_reads,
            outdir,
            long_reads=long_reads,
            long_read_type=read_type,
            threads=threads,
            memory_gb=memory_gb,
        )
        used_engine = "spades"
        notes.append(
            "assembled both read types together with SPAdes. If this comes out "
            "broken up, try mode='long_then_polish', which often does better on an "
            "AT rich genome read with nanopore"
        )
    elif mode == "long_then_polish":
        from .polish import polish

        contigs = assemble_with_flye(
            long_reads, outdir, genome_size, read_type, threads
        )
        used_engine = "flye"
        if short_reads:
            polished: list[Record] = []
            for contig in contigs:
                result = polish(
                    contig.seq,
                    short_reads,
                    read_type="illumina",
                    rounds=1,
                    threads=threads,
                    circular=False,
                    workdir=outdir / "polish" / contig.name,
                )
                polished.append(Record(contig.name, result.sequence))
            contigs = polished
            notes.append("polished the long read assembly with the short reads")
        else:
            notes.append("no short reads were given, so nothing was polished")
    else:
        raise ValueError(f"unknown mode {mode!r}")

    before = len(contigs)
    contigs = [c for c in contigs if len(c.seq) >= min_contig_bp]
    if before != len(contigs):
        notes.append(
            f"dropped {before - len(contigs)} contigs shorter than {min_contig_bp} bp"
        )

    contigs.sort(key=lambda c: len(c.seq), reverse=True)
    contigs = [Record(f"nuclear_contig_{i + 1}", c.seq) for i, c in enumerate(contigs)]

    return NuclearAssembly(
        contigs=contigs,
        engine=used_engine,
        mode=mode,
        stats=summarise(contigs),
        notes=notes,
    )


def available_engines() -> dict[str, bool]:
    """Which external assemblers are installed."""
    return {"flye": have("flye"), "spades": have("spades.py"), "nagasm": True}
