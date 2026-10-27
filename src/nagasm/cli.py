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


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    handlers = {
        "assemble": _cmd_assemble,
        "triage": _cmd_triage,
        "stats": _cmd_stats,
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
