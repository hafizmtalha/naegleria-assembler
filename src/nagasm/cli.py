"""The `nagasm` command.

Three subcommands:

* `assemble` runs the whole thing, from reads to an assembly and a report.
* `triage` only sorts the reads into compartments and draws the depth plot. It
  is fast, and it is the right first thing to run on a new dataset, because it
  tells you whether the compartments are really there before you spend hours on
  an assembly.
* `stats` describes an assembly that already exists.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from . import __version__


def _add_common(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "-o", "--outdir", default="results", help="where to write everything (default: results)"
    )
    parser.add_argument(
        "-t", "--threads", type=int, default=4, help="CPU threads to use (default: 4)"
    )
    parser.add_argument(
        "-k", type=int, default=21, help="k-mer length for sorting reads by depth (default: 21)"
    )
    parser.add_argument(
        "--scale",
        type=int,
        default=100,
        help="keep about one k-mer in this many when counting. Larger uses less "
        "memory and is slightly rougher (default: 100)",
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="nagasm",
        description=(
            "An assembly pipeline built around the way Naegleria keeps its DNA. "
            "It sorts reads into compartments by depth before assembling "
            "anything, so that the small circles held in thousands of copies do "
            "not wreck the assembly of the chromosomes. Chromosomes go to Flye or "
            "SPAdes. The sorting, the small circles and the read support check "
            "are its own."
        ),
    )
    parser.add_argument("--version", action="version", version=f"nagasm {__version__}")
    sub = parser.add_subparsers(dest="command", required=True)

    # --- assemble ---
    run = sub.add_parser("assemble", help="run the whole pipeline")
    run.add_argument(
        "-l",
        "--long",
        nargs="*",
        default=[],
        metavar="FILE",
        help="long read files, FASTA or FASTQ, gzipped is fine",
    )
    run.add_argument(
        "-s",
        "--short",
        nargs="*",
        default=[],
        metavar="FILE",
        help="short read files, FASTA or FASTQ, gzipped is fine",
    )
    _add_common(run)
    run.add_argument(
        "--read-type",
        default="ont",
        choices=["ont", "ont-hq", "hifi", "pacbio", "illumina"],
        help="what the long reads came off (default: ont)",
    )
    run.add_argument(
        "--mode",
        default=None,
        choices=["long", "short", "hybrid", "long_then_polish"],
        help="how to assemble the nuclear compartment. The default is worked out "
        "from which read files you gave",
    )
    run.add_argument(
        "--engine",
        default="auto",
        choices=["auto", "nagasm"],
        help="auto uses Flye or SPAdes for the nuclear assembly. nagasm uses the "
        "assembler in this package, which is here to be read rather than to "
        "compete (default: auto)",
    )
    run.add_argument("--memory", type=int, default=12, help="memory limit in GB for SPAdes")
    run.add_argument(
        "--min-contig", type=int, default=500, help="drop contigs shorter than this"
    )
    run.add_argument(
        "--polish-rounds", type=int, default=2, help="rounds of racon polishing (default: 2)"
    )
    run.add_argument("--skip-nuclear", action="store_true", help="only do the high copy circles")
    run.add_argument("--skip-polish", action="store_true", help="leave the drafts unpolished")

    pics = run.add_argument_group("pictures drawn when the run finishes")
    pics.add_argument(
        "--plots",
        default="auto",
        choices=["auto", "all", "none"],
        help="auto draws a circle per circular element, a karyotype and a "
        "contiguity curve. all adds coverage and a dotplot, at the cost of a "
        "minimap2 run each (default: auto)",
    )
    pics.add_argument(
        "--plot-theme",
        default="light",
        choices=["light", "dark", "print"],
        help="print is greyscale, for a journal (default: light)",
    )
    pics.add_argument(
        "--plot-format",
        default="png",
        choices=["png", "svg", "pdf"],
        help="svg and pdf stay sharp at any size (default: png)",
    )

    # --- triage ---
    tri = sub.add_parser(
        "triage",
        help="sort reads into compartments and draw the depth plot, without assembling",
    )
    tri.add_argument("reads", nargs="+", help="read files, FASTA or FASTQ")
    _add_common(tri)

    # --- stats ---
    st = sub.add_parser("stats", help="describe an assembly that already exists")
    st.add_argument("fasta", help="an assembly in FASTA")

    # --- plot ---
    # Optional, all of it. Nothing in the assembly depends on any of this.
    pl = sub.add_parser(
        "plot",
        help="draw pictures of an assembly (optional). Circles, karyotypes, "
        "dotplots, coverage and contiguity",
    )
    pl.add_argument("fasta", help="the assembly or the elements to draw, FASTA")
    pl.add_argument(
        "--kind",
        default="auto",
        choices=["auto", "all", "circle", "karyotype", "dotplot", "coverage", "nx"],
        help="what to draw. auto picks circles for circular elements and a "
        "karyotype for anything else (default: auto)",
    )
    pl.add_argument(
        "-o", "--outdir", default="figures", help="where to write them (default: figures)"
    )
    pl.add_argument(
        "-l",
        "--long",
        nargs="*",
        default=[],
        help="reads, for the depth ring, the coverage plot and depth shading. "
        "Reference free, and the honest option",
    )
    pl.add_argument(
        "-r",
        "--reference",
        help="a reference FASTA, for identity rings and the dotplot. Needs "
        "minimap2, and needs you to already have the answer",
    )

    look = pl.add_argument_group("how it looks")
    look.add_argument(
        "--theme",
        default="light",
        choices=["light", "dark", "print"],
        help="print is greyscale, for a journal (default: light)",
    )
    look.add_argument(
        "--format",
        dest="image_format",
        default="png",
        choices=["png", "svg", "pdf"],
        help="svg and pdf stay sharp at any size (default: png)",
    )
    look.add_argument("--dpi", type=int, default=150, help="for png (default: 150)")
    look.add_argument("--title", default="", help="your own title instead of the filename")
    look.add_argument("--no-legend", action="store_true", help="drop the legend")
    look.add_argument(
        "--figsize",
        help="width,height in inches, for example 12,8. The default fits the content",
    )

    circ = pl.add_argument_group("circles")
    circ.add_argument(
        "--rings",
        help="comma separated subset of: "
        "gc, skew, cumulative, depth, identity, missing. Default is all of them",
    )
    circ.add_argument(
        "-w",
        "--window",
        type=int,
        help="window size in bp. The default aims for a few hundred windows "
        "whatever the length",
    )
    circ.add_argument(
        "--only-circular",
        action="store_true",
        help="skip anything not marked circular",
    )
    circ.add_argument(
        "--min-length", type=int, default=100, help="skip anything shorter (default: 100)"
    )

    karyo = pl.add_argument_group("karyotype")
    karyo.add_argument(
        "--colour-by",
        default="gc",
        choices=["gc", "depth", "shape", "none"],
        help="what shades the bars. depth needs reads (default: gc)",
    )
    karyo.add_argument(
        "--max-contigs", type=int, default=40, help="draw the largest this many (default: 40)"
    )
    karyo.add_argument(
        "--sort-by", default="length", choices=["length", "name"], help="(default: length)"
    )
    karyo.add_argument(
        "--highlight", nargs="*", default=[], help="contig names to outline and bold"
    )

    other = pl.add_argument_group("dotplot, coverage and contiguity")
    other.add_argument(
        "--min-block",
        type=int,
        default=200,
        help="shortest alignment to draw in a dotplot (default: 200)",
    )
    other.add_argument(
        "--log-depth", action="store_true", help="log scale on the coverage plot"
    )
    other.add_argument(
        "--compare",
        nargs="*",
        default=[],
        help="more assemblies to put on the contiguity curve",
    )

    pl.add_argument(
        "-t", "--threads", type=int, default=4, help="CPU threads to use (default: 4)"
    )
    pl.add_argument(
        "--read-type",
        default="ont",
        choices=["ont", "pacbio", "hifi"],
        help="what kind of long reads (default: ont)",
    )

    # --- circleplot, kept as the older name for plot --kind circle ---
    cp = sub.add_parser("circleplot", help="the same as: plot --kind circle")
    cp.add_argument("fasta")
    cp.add_argument("-o", "--outdir", default="figures")
    cp.add_argument("-l", "--long", nargs="*", default=[])
    cp.add_argument("-r", "--reference")
    cp.add_argument("-w", "--window", type=int)
    cp.add_argument("-t", "--threads", type=int, default=4)
    cp.add_argument("--read-type", default="ont", choices=["ont", "pacbio", "hifi"])
    cp.add_argument("--theme", default="light", choices=["light", "dark", "print"])
    cp.add_argument("--rings")
    cp.add_argument("--format", dest="image_format", default="png", choices=["png", "svg", "pdf"])
    cp.add_argument("--dpi", type=int, default=150)
    cp.add_argument("--no-legend", action="store_true")
    cp.add_argument("--only-circular", action="store_true")
    cp.add_argument("--min-length", type=int, default=100)
    cp.add_argument("--title", default="")

    return parser


def _cmd_assemble(args) -> int:
    from .pipeline import Settings, run

    if not args.long and not args.short:
        print("Give at least one read file with --long or --short.", file=sys.stderr)
        return 2

    settings = Settings(
        outdir=Path(args.outdir),
        k=args.k,
        scale=args.scale,
        threads=args.threads,
        memory_gb=args.memory,
        read_type=args.read_type,
        engine=args.engine,
        mode=args.mode,
        min_contig_bp=args.min_contig,
        polish_rounds=args.polish_rounds,
        skip_nuclear=args.skip_nuclear,
        skip_polish=args.skip_polish,
        plots=args.plots,
        plot_theme=args.plot_theme,
        plot_format=args.plot_format,
    )
    run(args.long, args.short, settings)
    return 0


def _cmd_triage(args) -> int:
    from . import report
    from .pipeline import load_reads
    from .triage import triage

    reads = load_reads(args.reads)
    if not reads:
        print("No reads were loaded.", file=sys.stderr)
        return 1

    print(f"Loaded {len(reads):,} reads.")
    result = triage(reads, k=args.k, scale=args.scale)

    outdir = Path(args.outdir)
    (outdir / "figures").mkdir(parents=True, exist_ok=True)

    print()
    header = f"{'compartment':<14}{'depth':>9}{'vs nuclear':>12}{'span':>13}{'reads':>10}{'GC':>7}"
    print(header)
    print("-" * len(header))
    for comp in result.compartments:
        rel = "1x" if comp.label == "nuclear" else f"{comp.depth_relative_to_nuclear:.0f}x"
        print(
            f"{comp.label:<14}{comp.peak_depth:>8.0f}x{rel:>12}"
            f"{comp.span_bp:>12,d} {comp.n_reads:>9,d}{comp.median_gc * 100:>6.1f}%"
        )
        if comp.hint:
            print(f"{'':<14}{comp.hint}")

    if result.warnings:
        print()
        for warning in result.warnings:
            print(f"note: {warning}")

    path = report.plot_depth_spectrum(result, outdir / "figures" / "depth_spectrum.png")
    report.plot_reads_depth_vs_gc(result, outdir / "figures" / "reads_depth_vs_gc.png")
    print()
    print(f"Depth plot written to {path}. Look at it before trusting the table above.")
    return 0


def _cmd_stats(args) -> int:
    from .seqio import read_fasta
    from .stats import summarise

    records = list(read_fasta(args.fasta))
    if not records:
        print("That file has no sequences in it.", file=sys.stderr)
        return 1
    for line in summarise(records).summary_lines():
        print(line)
    return 0


def _has_circular(fasta: Path) -> bool:
    from .seqio import read_sequences

    return any("circular" in r.name.lower() for r in read_sequences(fasta))


def _cmd_plot(args) -> int:
    from . import circleplot, plots
    from .pipeline import load_reads

    fasta = Path(args.fasta)
    outdir = Path(args.outdir)
    reads = load_reads(args.long) if args.long else None
    if args.long and not reads:
        print("No reads were loaded, so anything needing them is skipped.", file=sys.stderr)
        reads = None
    reference = Path(args.reference) if getattr(args, "reference", None) else None
    rings = [r.strip() for r in args.rings.split(",")] if getattr(args, "rings", None) else None
    figsize = None
    if getattr(args, "figsize", None):
        try:
            w, h = (float(v) for v in args.figsize.split(","))
            figsize = (w, h)
        except ValueError:
            print("--figsize wants two numbers, like 12,8", file=sys.stderr)
            return 2

    kind = getattr(args, "kind", "circle")
    if kind == "auto":
        kind = "circle" if _has_circular(fasta) else "karyotype"
    kinds = ["circle", "karyotype", "dotplot", "coverage", "nx"] if kind == "all" else [kind]

    written: list[Path] = []
    for which in kinds:
        try:
            if which == "circle":
                written += circleplot.plot_all(
                    fasta,
                    outdir,
                    reads=reads,
                    reference=reference,
                    read_type=args.read_type,
                    threads=args.threads,
                    window=getattr(args, "window", None),
                    rings=rings,
                    theme_name=args.theme,
                    legend=not args.no_legend,
                    image_format=args.image_format,
                    dpi=args.dpi,
                    only_circular=getattr(args, "only_circular", False),
                    min_length=getattr(args, "min_length", 100),
                )
            elif which == "karyotype":
                written.append(
                    plots.karyotype(
                        fasta,
                        outdir / f"karyotype.{args.image_format}",
                        colour_by=args.colour_by,
                        reads=reads,
                        max_contigs=args.max_contigs,
                        min_length=args.min_length,
                        sort_by=args.sort_by,
                        highlight=args.highlight,
                        theme_name=args.theme,
                        title=args.title,
                        window=getattr(args, "window", None),
                        read_type=args.read_type,
                        threads=args.threads,
                        figsize=figsize,
                        dpi=args.dpi,
                    )
                )
            elif which == "dotplot":
                made = plots.dotplot(
                    fasta,
                    outdir / f"dotplot.{args.image_format}",
                    reference=reference,
                    min_block=args.min_block,
                    theme_name=args.theme,
                    title=args.title,
                    threads=args.threads,
                    dpi=args.dpi,
                    figsize=figsize or (8.0, 8.0),
                )
                if made is None:
                    print("dotplot: nothing aligned, or minimap2 is missing.", file=sys.stderr)
                else:
                    written.append(made)
            elif which == "coverage":
                if not reads:
                    print("coverage needs reads. Pass --long.", file=sys.stderr)
                    continue
                made = plots.coverage_plot(
                    fasta,
                    reads,
                    outdir / f"coverage.{args.image_format}",
                    window=getattr(args, "window", None),
                    theme_name=args.theme,
                    title=args.title,
                    read_type=args.read_type,
                    threads=args.threads,
                    dpi=args.dpi,
                    log_scale=args.log_depth,
                )
                if made is None:
                    print("coverage: nothing mapped, or minimap2 is missing.", file=sys.stderr)
                else:
                    written.append(made)
            elif which == "nx":
                compare = {Path(c).stem: Path(c) for c in getattr(args, "compare", [])}
                written.append(
                    plots.nx_curve(
                        fasta,
                        outdir / f"contiguity.{args.image_format}",
                        compare=compare or None,
                        theme_name=args.theme,
                        title=args.title,
                        figsize=figsize or (8.0, 5.0),
                        dpi=args.dpi,
                    )
                )
        except ValueError as exc:
            print(f"{which}: {exc}", file=sys.stderr)

    if not written:
        print("Nothing was drawn.", file=sys.stderr)
        return 1
    for path in written:
        print(f"wrote {path}")
    print(f"\n{len(written)} picture(s).")
    return 0


def _cmd_circleplot(args) -> int:
    args.kind = "circle"
    args.colour_by = "gc"
    args.max_contigs = 40
    args.sort_by = "length"
    args.highlight = []
    args.min_block = 200
    args.log_depth = False
    args.compare = []
    args.figsize = None
    return _cmd_plot(args)


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    handlers = {
        "assemble": _cmd_assemble,
        "triage": _cmd_triage,
        "stats": _cmd_stats,
        "plot": _cmd_plot,
        "circleplot": _cmd_circleplot,
    }
    try:
        return handlers[args.command](args)
    except KeyboardInterrupt:
        print("\nStopped.", file=sys.stderr)
        return 130
    except (ValueError, FileNotFoundError) as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
