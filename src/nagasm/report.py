"""Figures and a written report.

The depth plot is the one that matters. Everything this package does rests on
there being separate humps in the k-mer depth histogram, so the first thing
anyone should do with a new dataset is look at that plot and check the humps are
really there. If they are not, the compartment table is guesswork and the report
says so.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

# A colourblind friendly set, dark enough to read on white.
_COLOURS = ["#4477AA", "#EE6677", "#228833", "#CCBB44", "#66CCEE", "#AA3377", "#BBBBBB"]


def _style(ax) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(alpha=0.25, linewidth=0.6)
    ax.set_axisbelow(True)


def plot_depth_spectrum(result, path: Path, title: str = "") -> Path:
    """The k-mer depth histogram with the compartments marked.

    Both axes are logarithmic, and that is the point of the figure rather than a
    detail of it. On linear axes the nuclear hump fills the plot and the rDNA
    hump, which is a thousand times deeper and a thousand times smaller, is
    invisible. On log axes all the compartments are legible at once.
    """
    hist = np.asarray(result.histogram, dtype=float)
    depths = np.arange(hist.size)

    fig, ax = plt.subplots(figsize=(9, 5))
    keep = (depths >= 1) & (hist > 0)
    ax.plot(depths[keep], hist[keep], color="#333333", linewidth=1.0)

    ax.axvline(
        result.error_threshold,
        color="#999999",
        linestyle=":",
        linewidth=1.2,
        label=f"error cut ({result.error_threshold})",
    )

    for i, comp in enumerate(result.compartments):
        colour = _COLOURS[i % len(_COLOURS)]
        ax.axvspan(comp.low_depth, comp.high_depth, color=colour, alpha=0.13)
        ax.axvline(comp.peak_depth, color=colour, linewidth=1.6)
        label = f"{comp.label}\n{comp.peak_depth:.0f}x, {_bp(comp.span_bp)}"
        ax.annotate(
            label,
            xy=(comp.peak_depth, 1),
            xycoords=("data", "axes fraction"),
            xytext=(0, -14),
            textcoords="offset points",
            ha="center",
            va="top",
            fontsize=8,
            color=colour,
            fontweight="bold",
        )

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("How many times a k-mer is seen")
    ax.set_ylabel("How many distinct k-mers")
    ax.set_title(title or "Depth spectrum, with compartments marked")
    ax.legend(frameon=False, fontsize=8, loc="upper right")
    _style(ax)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return path


def plot_reads_depth_vs_gc(result, path: Path, max_points: int = 20_000) -> Path:
    """Every read as a dot, depth against GC, coloured by compartment.

    Depth alone decides the split, so this figure is a check rather than an
    input. Compartments that separate on GC as well are compartments you can
    trust. A cloud at a different GC from everything else is usually another
    organism in the culture rather than a part of this one.
    """
    depth = np.asarray(result.read_depth)
    gc = np.asarray(result.read_gc) * 100
    assignment = np.asarray(result.assignment)

    if depth.size > max_points:
        rng = np.random.default_rng(0)
        pick = rng.choice(depth.size, max_points, replace=False)
        depth, gc, assignment = depth[pick], gc[pick], assignment[pick]

    fig, ax = plt.subplots(figsize=(8, 5.5))
    unplaced = assignment < 0
    if unplaced.any():
        ax.scatter(
            np.maximum(depth[unplaced], 0.5), gc[unplaced], s=4, c="#DDDDDD",
            label="not placed", linewidths=0,
        )
    for i, comp in enumerate(result.compartments):
        sel = assignment == i
        if not sel.any():
            continue
        ax.scatter(
            np.maximum(depth[sel], 0.5),
            gc[sel],
            s=5,
            c=_COLOURS[i % len(_COLOURS)],
            label=f"{comp.label} ({comp.n_reads:,} reads)",
            linewidths=0,
            alpha=0.5,
        )

    ax.set_xscale("log")
    ax.set_xlabel("Read depth (median count of its k-mers)")
    ax.set_ylabel("GC of the read (%)")
    ax.set_title("Where each read sits")
    ax.legend(frameon=False, fontsize=8, markerscale=2)
    _style(ax)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return path


def plot_read_lengths(result, path: Path) -> Path:
    """Read length by compartment.

    Worth a look because a compartment made of unusually short reads is often an
    artefact rather than a real part of the genome.
    """
    assignment = np.asarray(result.assignment)
    lengths = np.asarray(result.read_length)

    groups = []
    labels = []
    for i, comp in enumerate(result.compartments):
        sel = assignment == i
        if sel.sum() >= 10:
            groups.append(lengths[sel])
            labels.append(f"{comp.label}\n{comp.n_reads:,}")
    if not groups:
        groups, labels = [lengths], ["all reads"]

    fig, ax = plt.subplots(figsize=(1.6 * len(groups) + 3, 5))
    parts = ax.boxplot(groups, labels=labels, patch_artist=True, showfliers=False)
    for patch, colour in zip(parts["boxes"], _COLOURS):
        patch.set_facecolor(colour)
        patch.set_alpha(0.55)
        patch.set_edgecolor("#333333")
    for whisker in parts["whiskers"] + parts["caps"]:
        whisker.set_color("#333333")
    for median in parts["medians"]:
        median.set_color("#111111")
        median.set_linewidth(1.6)

    ax.set_yscale("log")
    ax.set_ylabel("Read length (bp)")
    ax.set_title("Read length by compartment")
    _style(ax)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return path


def _bp(n: int) -> str:
    """A length written the way a person would say it."""
    if n >= 1_000_000:
        return f"{n / 1_000_000:.1f} Mb"
    if n >= 1_000:
        return f"{n / 1_000:.1f} kb"
    return f"{n} bp"


def write_summary(payload: dict[str, Any], path: Path) -> Path:
    """Everything the run worked out, as JSON, for scripting against."""
    path.write_text(json.dumps(payload, indent=2, default=str))
    return path


def write_report(payload: dict[str, Any], path: Path) -> Path:
    """A written summary of the run, in plain sentences."""
    lines: list[str] = []
    add = lines.append

    add("# Assembly report")
    add("")
    add(f"Run on {payload.get('started', 'an unrecorded date')}.")
    add("")

    reads = payload.get("reads", {})
    if reads:
        add("## The reads")
        add("")
        for name, info in reads.items():
            if not info.get("n_reads"):
                continue
            add(
                f"* {name}: {info['n_reads']:,} reads, {_bp(info['total_bp'])} in total, "
                f"read N50 {info['n50_bp']:,} bp, longest {info['longest_bp']:,} bp."
            )
        add("")

    triage = payload.get("triage")
    if triage:
        add("## Compartments")
        add("")
        add(
            "Reads were sorted by how often their k-mers turn up across the whole "
            "read set. Each group below is a hump in that distribution. The span is "
            "how much distinct sequence sits at that depth, which is what separates "
            "a large genome at ordinary depth from a small element held in many "
            "copies."
        )
        add("")
        add("| Compartment | Depth | Depth vs nuclear | Span | Reads | GC |")
        add("|---|---|---|---|---|---|")
        for comp in triage["compartments"]:
            rel = comp.get("depth_relative_to_nuclear")
            rel_text = "1x" if comp["label"] == "nuclear" else f"{rel:.0f}x"
            add(
                f"| {comp['label']} | {comp['peak_depth']:.0f}x | {rel_text} | "
                f"{_bp(comp['span_bp'])} | {comp['n_reads']:,} | "
                f"{comp['median_gc'] * 100:.1f}% |"
            )
        add("")
        for comp in triage["compartments"]:
            if comp.get("hint"):
                add(f"* **{comp['label']}**: {comp['hint']}.")
        add("")
        if triage.get("warnings"):
            add("### Things worth checking")
            add("")
            for warning in triage["warnings"]:
                add(f"* {warning}")
            add("")

    elements = payload.get("elements", [])
    if elements:
        add("## High copy elements")
        add("")
        for element in elements:
            shape = "a circle" if element["is_circular"] else "a line"
            add(
                f"### {element['name']}: {element['length']:,} bp, {shape}"
            )
            add("")
            if element.get("repeat_content") is not None:
                share = element["repeat_content"]
                add(
                    f"* About {share:.0%} of it sits inside sequence that repeats "
                    f"elsewhere in the same element."
                )
            for note in element.get("notes", []):
                add(f"* {note}.")
            if element.get("resolvability"):
                add(f"* {element['resolvability']}")
            add("")

    nuclear = payload.get("nuclear")
    if nuclear:
        add("## Nuclear assembly")
        add("")
        if nuclear.get("error"):
            # The nuclear step records an error instead of an engine when it
            # could not run at all, usually a missing Flye. Assuming the engine
            # key was always there crashed the whole report writer, which threw
            # away a run that had already found and polished every element.
            add(f"This step did not finish: {nuclear['error']}")
            add("")
            add(
                "The compartments above were still sorted and the high copy "
                "elements above were still assembled, so the rest of this report "
                "stands. Only the chromosomes are missing."
            )
            add("")
        else:
            add(
                f"Assembled with {nuclear.get('engine', 'an unrecorded engine')} in "
                f"{nuclear.get('mode', 'an unrecorded')} mode."
            )
            add("")
        stats = nuclear.get("stats", {})
        if stats:
            add(f"* Contigs: {stats['n_contigs']:,}")
            add(f"* Total length: {stats['total_bp']:,} bp")
            add(f"* Longest contig: {stats['longest_bp']:,} bp")
            add(f"* N50: {stats['n50']:,} bp (L50 {stats['l50']:,})")
            add(f"* GC: {stats['gc_percent']:.1f}%")
        add("")
        for note in nuclear.get("notes", []):
            add(f"* {note}.")
        add("")

    busco = payload.get("completeness")
    if busco:
        add("## Gene content")
        add("")
        if busco.get("error"):
            add(
                f"The completeness check was asked for but did not finish: "
                f"{busco['error']}"
            )
            add("")
        else:
            version = busco.get("busco_version")
            tool = f"BUSCO {version}" if version else "BUSCO"
            add(
                f"Scored with {tool} against "
                f"**{busco.get('lineage', 'an unnamed lineage')}**, "
                f"{busco.get('n_markers', 0):,} markers, on "
                f"`{busco.get('scored', 'the assembly')}`."
            )
            add("")
            add("| | Share | Count |")
            add("|---|---|---|")
            rows = (
                ("Complete", "complete_pct", "complete_n"),
                ("of those, single copy", "single_pct", "single_n"),
                ("of those, duplicated", "duplicated_pct", "duplicated_n"),
                ("Fragmented", "fragmented_pct", "fragmented_n"),
                ("Missing", "missing_pct", "missing_n"),
            )
            for label, pct_key, count_key in rows:
                add(
                    f"| {label} | {busco.get(pct_key, 0.0):.1f}% | "
                    f"{busco.get(count_key, 0):,} |"
                )
            add("")
            if busco.get("one_line"):
                add(f"BUSCO's own one line summary: `{busco['one_line']}`")
                add("")
            for note in busco.get("notes", []):
                add(f"* {note}.")
            add("")

    add("## A note on reading these numbers")
    add("")
    add(
        "N50 rewards long contigs and says nothing about whether they are right. "
        "An assembly that glued two chromosomes together would score better than "
        "one that left them apart. Read the contig numbers next to the compartment "
        "table, not instead of it."
    )
    add("")

    path.write_text("\n".join(lines))
    return path
