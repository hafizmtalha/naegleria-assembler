"""The whole run, from reads to an assembly and a report.

The order of the steps is the argument this package is making, so it is worth
stating plainly:

1. Read the reads.
2. Sort them into compartments by depth.
3. Assemble each high copy compartment on its own and close it into a circle.
4. Assemble the nuclear compartment from what is left.
5. Polish everything.
6. Write the figures and the report.

Step 2 before step 4 is the whole idea. Doing it the other way round, or not at
all, is what makes a general assembler struggle on this organism.
"""

from __future__ import annotations

import datetime as _dt
from dataclasses import dataclass, field
from pathlib import Path

from . import circular, nuclear as nuclear_mod, polish as polish_mod, report, stats, triage
from .seqio import Record, read_sequences, write_fasta


@dataclass
class Settings:
    """Everything the run can be told to do differently."""

    outdir: Path = Path("results")
    k: int = 21
    scale: int = 100
    threads: int = 4
    memory_gb: int = 12
    read_type: str = "ont"
    engine: str = "auto"
    mode: str | None = None
    min_contig_bp: int = 500
    polish_rounds: int = 2
    circle_depth: float = 60.0
    circle_k: int = 15
    max_elements: int = 4
    # Reads this many times deeper than the nuclear peak are assembled on their
    # own even when they never formed a hump. See
    # TriageResult.deep_outlier_indices for why this is needed.
    deep_outlier_factor: float = 4.0
    # How many circles to report from the deep outlier search. That search is a
    # heuristic rather than a measurement, and most of what it turns up beyond
    # the first hit is nuclear repeat that the walk happened to close. Two keeps
    # the real find without burying it.
    max_deep_elements: int = 2
    skip_nuclear: bool = False
    skip_polish: bool = False


@dataclass
class RunResult:
    settings: Settings
    payload: dict = field(default_factory=dict)
    outputs: dict[str, Path] = field(default_factory=dict)


def load_reads(paths: list[str | Path]) -> list[str]:
    """Read every sequence from every file given."""
    out: list[str] = []
    for path in paths:
        out.extend(record.seq for record in read_sequences(path))
    return out


def run(
    long_read_paths: list[str | Path] | None = None,
    short_read_paths: list[str | Path] | None = None,
    settings: Settings | None = None,
    log=print,
) -> RunResult:
    """Run the whole thing."""
    settings = settings or Settings()
    outdir = Path(settings.outdir)
    (outdir / "figures").mkdir(parents=True, exist_ok=True)
    (outdir / "assembly").mkdir(parents=True, exist_ok=True)

    long_reads = load_reads(long_read_paths or [])
    short_reads = load_reads(short_read_paths or [])
    if not long_reads and not short_reads:
        raise ValueError("no reads were loaded, so there is nothing to do")

    payload: dict = {
        "started": _dt.datetime.now().strftime("%Y-%m-%d %H:%M"),
        "settings": {
            key: str(value) for key, value in settings.__dict__.items()
        },
        "reads": {
            "long": stats.read_length_summary(long_reads),
            "short": stats.read_length_summary(short_reads),
        },
    }
    outputs: dict[str, Path] = {}

    log(
        f"Loaded {len(long_reads):,} long reads and {len(short_reads):,} short reads."
    )

    # --- step 2: sort the reads into compartments ------------------------
    #
    # Triage runs on the long reads when there are any. They span more of each
    # element, so a long read's depth describes where it came from far better
    # than a short read's does.
    triage_reads = long_reads or short_reads
    log("Sorting reads into compartments by k-mer depth.")
    triage_result = triage.triage(triage_reads, k=settings.k, scale=settings.scale)

    payload["triage"] = {
        "k": triage_result.k,
        "scale": triage_result.scale,
        "hpc": triage_result.hpc,
        "error_threshold": triage_result.error_threshold,
        "compartments": [c.as_dict() for c in triage_result.compartments],
        "warnings": triage_result.warnings,
    }
    for comp in triage_result.compartments:
        log(
            f"  {comp.label}: {comp.peak_depth:.0f}x deep, about "
            f"{comp.span_bp:,} bp of distinct sequence, {comp.n_reads:,} reads"
        )
    for warning in triage_result.warnings:
        log(f"  note: {warning}")

    outputs["depth_spectrum"] = report.plot_depth_spectrum(
        triage_result, outdir / "figures" / "depth_spectrum.png"
    )
    outputs["reads_depth_vs_gc"] = report.plot_reads_depth_vs_gc(
        triage_result, outdir / "figures" / "reads_depth_vs_gc.png"
    )
    outputs["read_lengths"] = report.plot_read_lengths(
        triage_result, outdir / "figures" / "read_lengths.png"
    )

    # --- step 3: the high copy circles -----------------------------------
    elements: list[dict] = []
    element_records: list[Record] = []
    high_copy = sorted(triage_result.high_copy, key=lambda c: c.peak_depth, reverse=True)
    for comp in high_copy[: settings.max_elements]:
        idx = triage_result.read_indices(comp.label)
        comp_reads = [triage_reads[i] for i in idx]
        if len(comp_reads) < 5:
            log(f"  {comp.label}: only {len(comp_reads)} reads, skipping")
            continue

        log(f"Assembling {comp.label} ({len(comp_reads):,} reads).")
        # A compartment can hold more than one molecule, so keep going until
        # nothing more comes out rather than stopping at the first element.
        assemblies = circular.assemble_elements(
            comp_reads,
            expected_bp=comp.span_bp,
            max_elements=settings.max_elements,
            target_depth=settings.circle_depth,
            k=settings.circle_k,
        )
        if not assemblies:
            log(f"  {comp.label}: nothing assembled")
            continue

        for part, asm in enumerate(assemblies, start=1):
            sequence = asm.sequence
            if not settings.skip_polish and polish_mod.have("minimap2") and polish_mod.have("racon"):
                log(f"  polishing {comp.label} part {part} ({len(sequence):,} bp).")
                polished = polish_mod.polish(
                    sequence,
                    comp_reads,
                    read_type=settings.read_type,
                    rounds=settings.polish_rounds,
                    threads=settings.threads,
                    circular=asm.is_circular,
                    workdir=outdir / "work" / f"polish_{comp.label}_{part}",
                )
                sequence = polished.sequence
                asm.notes.extend(polished.notes)
                if asm.is_circular:
                    sequence = circular.rotate_to_canonical_start(sequence)

            suffix = "circular" if asm.is_circular else "linear"
            name = f"{comp.label}_{part}_{suffix}"
            element_records.append(Record(name, sequence))
            elements.append(
                {
                    "name": name,
                    "compartment": comp.label,
                    "length": len(sequence),
                    "is_circular": asm.is_circular,
                    "mean_coverage": round(asm.mean_coverage, 1),
                    "reads_used": asm.n_reads_used,
                    "notes": asm.notes,
                    "repeat_content": round(circular.repeat_content(sequence, 50), 3),
                    "resolvability": circular.resolvability_note(sequence, comp_reads),
                }
            )
            log(
                f"  {name}: {len(sequence):,} bp, "
                f"{'closed into a circle' if asm.is_circular else 'left as a line'}"
            )

    # Elements hiding in the tail above the nuclear peak, with no hump of their
    # own. On Naegleria this is where the rDNA circle usually is.
    outlier_idx = triage_result.deep_outlier_indices(settings.deep_outlier_factor)
    if len(outlier_idx) >= 30:
        outlier_reads = [triage_reads[i] for i in outlier_idx]
        log(
            f"Looking for elements among the {len(outlier_reads):,} reads more than "
            f"{settings.deep_outlier_factor:.0f} times deeper than the nuclear peak."
        )
        for n, asm in enumerate(
            circular.assemble_elements(
                outlier_reads,
                expected_bp=0,
                max_elements=settings.max_deep_elements,
                target_depth=settings.circle_depth,
                k=settings.circle_k,
            ),
            start=1,
        ):
            if not asm.is_circular:
                # Nuclear repeats make up most of this set and they come out as
                # lines. Only the circles are worth reporting as elements.
                continue
            sequence = asm.sequence
            if not settings.skip_polish and polish_mod.have("minimap2") and polish_mod.have("racon"):
                polished = polish_mod.polish(
                    sequence,
                    outlier_reads,
                    read_type=settings.read_type,
                    rounds=settings.polish_rounds,
                    threads=settings.threads,
                    circular=True,
                    workdir=outdir / "work" / f"polish_deep_{n}",
                )
                sequence = circular.rotate_to_canonical_start(polished.sequence)
                asm.notes.extend(polished.notes)

            # The deep reads include the high copy compartments, so an element
            # already assembled from one of those turns up here again.
            #
            # This comparison has to happen after polishing, not before. An
            # unpolished draft carries the error rate of the reads, so its
            # anchor k-mers barely overlap those of a polished sequence even
            # when both describe the same molecule. Comparing the two states
            # let the mitochondrial genome through a second time.
            if any(circular.same_element(sequence, rec.seq) for rec in element_records):
                log("  this is an element already assembled from a compartment, skipping")
                continue

            name = f"deep_element_{n}_circular"
            element_records.append(Record(name, sequence))
            elements.append(
                {
                    "name": name,
                    "compartment": "above the nuclear peak, no hump of its own",
                    "length": len(sequence),
                    "is_circular": True,
                    "mean_coverage": round(asm.mean_coverage, 1),
                    "reads_used": asm.n_reads_used,
                    "notes": asm.notes,
                    "repeat_content": round(circular.repeat_content(sequence, 50), 3),
                    "resolvability": circular.resolvability_note(sequence, outlier_reads),
                }
            )
            log(f"  {name}: {len(sequence):,} bp, closed into a circle")

    payload["elements"] = elements
    if element_records:
        path = outdir / "assembly" / "high_copy_elements.fasta"
        write_fasta(path, element_records)
        outputs["elements"] = path

    # --- step 4: the nuclear compartment ---------------------------------
    if settings.skip_nuclear:
        log("Skipping the nuclear assembly, as asked.")
    else:
        nuc = triage_result.nuclear
        if nuc is None:
            log("No nuclear compartment was found, so there is nothing to assemble.")
        else:
            nuclear_long = _reads_for(long_reads, triage_result, nuc.label) if long_reads else []
            # Short reads were not triaged, so they all go in. The high copy
            # elements are already assembled separately, and leaving their short
            # reads in costs a little time rather than correctness.
            nuclear_short = short_reads

            log(
                f"Assembling the nuclear compartment "
                f"({len(nuclear_long):,} long, {len(nuclear_short):,} short reads, "
                f"estimated {nuc.span_bp:,} bp)."
            )
            try:
                result = nuclear_mod.assemble_nuclear(
                    long_reads=nuclear_long,
                    short_reads=nuclear_short,
                    outdir=outdir / "work" / "nuclear",
                    genome_size=nuc.span_bp,
                    read_type=settings.read_type,
                    mode=settings.mode,
                    engine=settings.engine,
                    threads=settings.threads,
                    memory_gb=settings.memory_gb,
                    min_contig_bp=settings.min_contig_bp,
                )
            except (RuntimeError, polish_mod.ToolMissing) as exc:
                log(f"  the nuclear assembly did not finish: {exc}")
                payload["nuclear"] = {"error": str(exc)}
            else:
                path = outdir / "assembly" / "nuclear.fasta"
                write_fasta(path, result.contigs)
                outputs["nuclear"] = path
                payload["nuclear"] = {
                    "engine": result.engine,
                    "mode": result.mode,
                    "stats": result.stats.as_dict() if result.stats else {},
                    "notes": result.notes,
                }
                if result.stats:
                    log(
                        f"  {result.stats.n_contigs:,} contigs, "
                        f"{result.stats.total_bp:,} bp, N50 {result.stats.n50:,} bp"
                    )

    # --- step 6: write everything out ------------------------------------
    outputs["summary"] = report.write_summary(payload, outdir / "summary.json")
    outputs["report"] = report.write_report(payload, outdir / "report.md")
    log(f"Report written to {outputs['report']}.")

    return RunResult(settings=settings, payload=payload, outputs=outputs)


def _reads_for(reads: list[str], result: triage.TriageResult, label: str) -> list[str]:
    """The reads assigned to one compartment, when triage ran on this read set."""
    if len(reads) != result.assignment.size:
        # Triage ran on a different set, so we cannot map indices across.
        return reads
    return [reads[i] for i in result.read_indices(label)]
