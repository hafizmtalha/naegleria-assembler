"""The whole run, from reads to an assembly and a report.

The order of the steps is the argument this package is making, so it is worth
stating plainly:

1. Read the reads.
2. Sort them into compartments by depth.
3. Assemble each high copy compartment on its own and close it into a circle.
4. Assemble the nuclear compartment from what is left.
5. Polish everything.
6. Score the nuclear contigs for gene content, if asked to.
7. Write the figures and the report.

Step 2 before step 4 is the whole idea. Doing it the other way round, or not at
all, is what makes a general assembler struggle on this organism.
"""

from __future__ import annotations

import datetime as _dt
from dataclasses import dataclass, field
from pathlib import Path

from . import (
    circular,
    nuclear as nuclear_mod,
    polish as polish_mod,
    report,
    stats,
    support,
    triage,
)
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
    # What to draw when the run finishes. "auto" draws a circle for every
    # circular element plus a karyotype and a contiguity curve, which is what
    # most people want and nobody remembers to ask for. "none" draws nothing.
    # "all" adds the coverage plot and the dotplot, which cost a minimap2 run
    # each.
    plots: str = "auto"
    plot_theme: str = "light"
    plot_format: str = "png"
    # Which BUSCO lineage to score the nuclear contigs against, or None to skip
    # the check entirely, which is the default. It is off because BUSCO is not in
    # environment.yml and because it downloads its lineage the first time, and
    # nothing else in this pipeline needs the network. See completeness.py.
    busco: str | None = None
    busco_download_path: str | None = None


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

            # Ask whether the reads actually back this up before calling it an
            # element. A deep compartment in a bacterium is usually a repeat
            # inside the chromosome rather than a molecule of its own, and the
            # assembler will happily close a circle out of it.
            # Measured against every read, not just the reads in this
            # compartment. The compartment was chosen for being deep, which
            # means chosen for reads that are mostly element sequence, so
            # asking it to judge its own element is a rigged jury. The reads
            # that give a chromosomal repeat away are exactly the ones that
            # straddle its edges, and those sit at ordinary depth and were
            # left out of the compartment. On real E. coli data this was the
            # difference between a repeat scoring 0.23 and scoring 0.10.
            element_support = support.measure(
                sequence,
                triage_reads,
                is_circular=asm.is_circular,
                read_type=settings.read_type,
                threads=settings.threads,
                workdir=outdir / "work" / f"support_{comp.label}_{part}",
            )
            supported, why = support.verdict(element_support)
            asm.notes.append(why)
            if not supported:
                log(f"  {name}: not supported by the reads, {why}")

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
                    "supported": supported,
                    "support": element_support.as_dict(),
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

            # These reads were collected because they are deep, not because
            # they form a hump of their own, so this path is the one most
            # likely to hand back a repeat dressed up as a circle.
            element_support = support.measure(
                sequence,
                triage_reads,
                is_circular=True,
                read_type=settings.read_type,
                threads=settings.threads,
                workdir=outdir / "work" / f"support_deep_{n}",
            )
            supported, why = support.verdict(element_support)
            asm.notes.append(why)
            if not supported:
                log(f"  {name}: not supported by the reads, {why}")

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
                    "supported": supported,
                    "support": element_support.as_dict(),
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

            # Take out any read that matches an element we already assembled.
            #
            # Sorting by depth gets most of them, but it does not have to miss
            # many to matter. On the real data 97 percent of the mitochondrial
            # reads went to the right compartment, and the 3 percent left behind
            # still amounted to sixty fold coverage of a 49 kb molecule, which
            # was plenty for the nuclear assembler to build the whole thing a
            # second time. A 49 kb mitochondrial contig sat in the middle of the
            # chromosomes, which is exactly what this package exists to avoid.
            #
            # Now that the elements are assembled, we can do better than a depth
            # cutoff and ask directly whether a read looks like one of them.
            # Only elements the reads actually support are taken out. An
            # element that failed that check is most likely a repeat living
            # inside the chromosome, and pulling its reads out would strip
            # real sequence from the assembly it belongs to.
            supported_names = {
                e["name"] for e in elements if e.get("supported", True)
            }
            keep_out = [r for r in element_records if r.name in supported_names]
            if keep_out and nuclear_long:
                before = len(nuclear_long)
                for record in keep_out:
                    if not nuclear_long:
                        break
                    belongs = circular.reads_from_element(
                        nuclear_long, record.seq, is_circular="circular" in record.name
                    )
                    nuclear_long = [
                        r for r, hit in zip(nuclear_long, belongs) if not hit
                    ]
                removed = before - len(nuclear_long)
                if removed:
                    log(
                        f"  removed {removed:,} reads that match an element already "
                        f"assembled, leaving {len(nuclear_long):,}"
                    )
                    payload.setdefault("nuclear_prefilter", {})[
                        "reads_removed_matching_elements"
                    ] = removed
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

    # --- step 6: optional completeness check -----------------------------
    _score_completeness(settings, outdir, outputs, payload, log)

    # --- step 7: write everything out ------------------------------------
    outputs["summary"] = report.write_summary(payload, outdir / "summary.json")
    outputs["report"] = report.write_report(payload, outdir / "report.md")
    log(f"Report written to {outputs['report']}.")

    drawn = _draw_plots(settings, outdir, outputs, payload, triage_reads, log)
    if drawn:
        outputs["plots"] = drawn

    return RunResult(settings=settings, payload=payload, outputs=outputs)


def _score_completeness(settings, outdir: Path, outputs: dict, payload: dict, log) -> None:
    """Run BUSCO on the nuclear contigs, when a lineage was asked for.

    Off unless `--busco <lineage>` is given. Like the figures, it never fails the
    run: a missing BUSCO, a lineage that will not download, or a crash inside it
    all end as a note in the report. Losing a finished assembly to an optional
    quality check would be absurd.

    It scores `nuclear.fasta` and nothing else, on purpose. Pointing BUSCO at the
    high copy elements as well would report duplicated markers that are not
    duplications, they are two kinds of DNA in one file. completeness.py explains
    this at more length.
    """
    if not settings.busco:
        return

    nuclear = outputs.get("nuclear")
    if not nuclear or not Path(nuclear).exists():
        log("  no nuclear contigs to score, so BUSCO was skipped")
        payload["completeness"] = {"error": "there were no nuclear contigs to score"}
        return

    from . import completeness as completeness_mod

    log(f"Scoring the nuclear contigs against {settings.busco}. BUSCO takes a while.")
    try:
        result = completeness_mod.run_busco(
            Path(nuclear),
            outdir / "completeness",
            lineage=settings.busco,
            threads=settings.threads,
            download_path=settings.busco_download_path,
        )
    except Exception as exc:
        log(f"  the completeness check did not finish: {exc}")
        payload["completeness"] = {"error": str(exc), "lineage": settings.busco}
        return

    payload["completeness"] = result.as_dict()
    outputs["completeness"] = Path(result.outdir)
    for line in result.summary_lines():
        log(f"  {line}")


def _draw_plots(settings, outdir: Path, outputs: dict, payload: dict, reads, log) -> list[Path]:
    """Draw the obvious pictures once the assembly is done.

    This runs by default, because the pictures are what most people want and
    nobody remembers to ask for them. It never fails the run. A missing
    matplotlib, a missing minimap2 or a sequence too short to draw all end in a
    note rather than an exception, since losing a finished assembly to a failed
    figure would be absurd.
    """
    if settings.plots == "none":
        return []

    try:
        from . import circleplot, plots as plots_mod
    except Exception as exc:  # matplotlib missing, most likely
        log(f"  no plots drawn: {exc}")
        return []

    figures = outdir / "figures"
    elements = outputs.get("elements")
    nuclear = outputs.get("nuclear")
    everything = settings.plots == "all"
    drawn: list[Path] = []

    # An element that failed the support check is still drawn, because it is
    # still reported, but the picture has to say so on its face.
    notes = {}
    for element in payload.get("elements", []):
        if element.get("supported", True):
            continue
        contained = element.get("support", {}).get("contained")
        warning = "NOT SUPPORTED BY THE READS"
        if contained is not None:
            warning += f", only {contained:.0%} of each read lies inside it"
        notes[element["name"]] = warning

    try:
        if elements and Path(elements).exists():
            drawn += circleplot.plot_all(
                Path(elements),
                figures,
                notes=notes,
                reads=reads if everything else None,
                read_type=settings.read_type,
                threads=settings.threads,
                theme_name=settings.plot_theme,
                image_format=settings.plot_format,
                only_circular=True,
            )
        if nuclear and Path(nuclear).exists():
            drawn.append(
                plots_mod.karyotype(
                    Path(nuclear),
                    figures / f"karyotype.{settings.plot_format}",
                    colour_by="gc",
                    theme_name=settings.plot_theme,
                    threads=settings.threads,
                )
            )
            drawn.append(
                plots_mod.nx_curve(
                    Path(nuclear),
                    figures / f"contiguity.{settings.plot_format}",
                    theme_name=settings.plot_theme,
                )
            )
        if everything and nuclear and Path(nuclear).exists() and reads:
            made = plots_mod.coverage_plot(
                Path(nuclear),
                reads,
                figures / f"coverage.{settings.plot_format}",
                theme_name=settings.plot_theme,
                read_type=settings.read_type,
                threads=settings.threads,
            )
            if made:
                drawn.append(made)
        if everything and elements and Path(elements).exists():
            made = plots_mod.dotplot(
                Path(elements),
                figures / f"dotplot.{settings.plot_format}",
                theme_name=settings.plot_theme,
                threads=settings.threads,
            )
            if made:
                drawn.append(made)
    except Exception as exc:
        log(f"  plots stopped early: {exc}")

    if drawn:
        log(f"Drew {len(drawn)} figure(s) in {figures}.")
        log("  Want them different, or want more of them?")
        log(f"    nagasm plot {elements or nuclear} --kind all -o {figures} --theme dark")
        log("  nagasm plot --help lists every option.")
    return drawn


def _reads_for(reads: list[str], result: triage.TriageResult, label: str) -> list[str]:
    """The reads assigned to one compartment, when triage ran on this read set."""
    if len(reads) != result.assignment.size:
        # Triage ran on a different set, so we cannot map indices across.
        return reads
    return [reads[i] for i in result.read_indices(label)]
