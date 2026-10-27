"""Pictures of an assembly, beyond the circles.

`circleplot.py` draws one circular element as a circle. This draws everything
else, and the two share a theme so a report does not look assembled from two
different tools.

Four kinds here:

* **karyotype**, every contig as a bar, scaled against each other, shaded by GC
  or by read depth. This is the picture that tells you at a glance whether you
  have a handful of chromosomes or a heap of rubble.
* **dotplot**, a sequence against itself or against a reference. Against itself
  it shows repeats, which is the thing that breaks assemblies, and the off
  diagonal lines are the repeat copies. Against a reference it shows
  rearrangement.
* **coverage**, read depth along a contig as a straight line rather than a ring.
  Better than the circular version when you want to read the numbers.
* **nx**, the contiguity curve. N50 is one point on this curve, and the curve
  tells you much more than the point does.

Everything here is optional. Nothing in the assembly imports it.
"""

from __future__ import annotations

import math
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path

import numpy as np

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap  # noqa: E402
from matplotlib.patches import Rectangle  # noqa: E402

from .seqio import Record, read_sequences, write_fasta, write_fastq
from .polish import have, _preset_for


# A colourblind friendly set, dark enough to read on white and light enough to
# read on near black.
PALETTE = ["#4477AA", "#EE6677", "#228833", "#CCBB44", "#66CCEE", "#AA3377", "#BBBBBB"]

THEMES = {
    "light": {
        "bg": "#FFFFFF",
        "fg": "#222222",
        "muted": "#666666",
        "grid": "#DDDDDD",
        "palette": PALETTE,
    },
    "dark": {
        "bg": "#14161A",
        "fg": "#E8E8E8",
        "muted": "#9AA0A6",
        "grid": "#2E3238",
        "palette": [
            "#6699DD",
            "#FF8899",
            "#55BB77",
            "#DDCC66",
            "#88DDFF",
            "#CC77AA",
            "#999999",
        ],
    },
    "print": {
        # Greys only, for a journal that will print it in black and white.
        "bg": "#FFFFFF",
        "fg": "#000000",
        "muted": "#555555",
        "grid": "#CCCCCC",
        "palette": ["#111111", "#555555", "#888888", "#AAAAAA", "#333333", "#777777", "#999999"],
    },
}

GRADIENTS = {
    "identity": ["#EE6677", "#CCBB44", "#228833"],
    "gc": ["#4477AA", "#EEEEEE", "#EE6677"],
    "depth": ["#F2F2F2", "#66CCEE", "#AA3377"],
    "mono": ["#EEEEEE", "#111111"],
}


def gradient(name: str) -> LinearSegmentedColormap:
    return LinearSegmentedColormap.from_list(name, GRADIENTS.get(name, GRADIENTS["mono"]))


def theme(name: str) -> dict:
    return THEMES.get(name, THEMES["light"])


def apply_theme(fig, ax_or_axes, name: str) -> dict:
    """Paint the figure and axes in a theme, and hand the theme back."""
    t = theme(name)
    fig.patch.set_facecolor(t["bg"])
    axes = ax_or_axes if isinstance(ax_or_axes, (list, tuple, np.ndarray)) else [ax_or_axes]
    for ax in np.ravel(np.asarray(axes, dtype=object)):
        if ax is None:
            continue
        ax.set_facecolor(t["bg"])
        for spine in ax.spines.values():
            spine.set_color(t["grid"])
        ax.tick_params(colors=t["muted"], labelcolor=t["muted"])
        ax.xaxis.label.set_color(t["fg"])
        ax.yaxis.label.set_color(t["fg"])
        ax.title.set_color(t["fg"])
    return t


def bp_label(n: float, total: float = 0) -> str:
    """A length label with enough precision that two of them never read alike."""
    scale = total or n
    if scale >= 1_000_000:
        return f"{n / 1_000_000:.2f} Mb" if scale < 10_000_000 else f"{n / 1_000_000:.1f} Mb"
    if scale >= 10_000:
        return f"{n / 1_000:.0f} kb"
    if scale >= 1_000:
        return f"{n / 1_000:.1f} kb"
    return f"{n:.0f} bp"


def safe_name(name: str) -> str:
    return "".join(c if c.isalnum() or c in "._-" else "_" for c in name)


def _windows(length: int, window: int) -> list[tuple[int, int]]:
    if window <= 0:
        raise ValueError("window must be positive")
    return [(s, min(s + window, length)) for s in range(0, length, window)]


def _gc(chunk: str) -> float:
    chunk = chunk.upper()
    gc = chunk.count("G") + chunk.count("C")
    at = chunk.count("A") + chunk.count("T")
    return gc / (gc + at) if gc + at else 0.0


# ---------------------------------------------------------------- karyotype


@dataclass
class Contig:
    name: str
    length: int
    sequence: str = ""
    is_circular: bool = False
    depth: float | None = None


def _load_contigs(fasta: Path, min_length: int = 0, keep_sequence: bool = True) -> list[Contig]:
    out = []
    for record in read_sequences(fasta):
        if len(record.seq) < min_length:
            continue
        out.append(
            Contig(
                name=record.name,
                length=len(record.seq),
                sequence=record.seq if keep_sequence else "",
                is_circular="circular" in record.name.lower(),
            )
        )
    return out


def karyotype(
    fasta: Path,
    path: Path,
    colour_by: str = "gc",
    reads: list[str] | None = None,
    max_contigs: int = 40,
    min_length: int = 0,
    sort_by: str = "length",
    highlight: list[str] | None = None,
    theme_name: str = "light",
    title: str = "",
    window: int | None = None,
    read_type: str = "ont",
    threads: int = 4,
    figsize: tuple[float, float] | None = None,
    dpi: int = 150,
) -> Path:
    """Every contig as a bar, drawn to the same scale.

    `colour_by` is "gc", "depth", "none" or "shape". Depth needs reads and
    minimap2, and quietly falls back to GC if either is missing, because a
    picture with one ring missing beats no picture at all.

    Only the largest `max_contigs` are drawn. A 65 contig assembly is readable.
    A 4,000 contig one is a black smear, and the honest thing is to draw the
    big ones and say in the title how many were left out.
    """
    contigs = _load_contigs(Path(fasta), min_length=min_length)
    if not contigs:
        raise ValueError("no contigs long enough to draw")

    total_all = sum(c.length for c in contigs)
    n_all = len(contigs)

    if sort_by == "name":
        contigs.sort(key=lambda c: c.name)
    else:
        contigs.sort(key=lambda c: c.length, reverse=True)
    hidden = max(0, len(contigs) - max_contigs)
    if sort_by == "name":
        biggest = set(
            c.name for c in sorted(contigs, key=lambda c: c.length, reverse=True)[:max_contigs]
        )
        contigs = [c for c in contigs if c.name in biggest]
    else:
        contigs = contigs[:max_contigs]

    longest = max(c.length for c in contigs)
    if window is None:
        window = max(50, longest // 500)

    if colour_by == "depth":
        depths = _depth_per_contig(
            contigs, reads, window, read_type=read_type, threads=threads
        )
        if depths is None:
            colour_by = "gc"
    else:
        depths = None

    # Fit the GC scale to the assembly. A fixed 25 to 75 percent band wastes
    # almost all of its range on an AT rich genome like Naegleria, where every
    # window sits near 35 percent and the whole picture comes out one colour.
    gc_range = (0.25, 0.75)
    if colour_by == "gc":
        sample = np.concatenate(
            [
                np.asarray([_gc(c.sequence[s:e]) for s, e in _windows(c.length, window)])
                for c in contigs[:10]
            ]
        )
        lo = float(np.percentile(sample, 2))
        hi = float(np.percentile(sample, 98))
        if hi - lo < 0.06:  # a floor, so noise is not magnified into stripes
            mid = (lo + hi) / 2
            lo, hi = mid - 0.03, mid + 0.03
        gc_range = (lo, hi)

    height = figsize[1] if figsize else max(3.0, 0.34 * len(contigs) + 2.0)
    width = figsize[0] if figsize else 11.0
    fig, ax = plt.subplots(figsize=(width, height))
    t = apply_theme(fig, ax, theme_name)
    highlight = set(highlight or [])

    cmap = gradient("gc" if colour_by == "gc" else "depth")
    bar_h = 0.52

    for i, contig in enumerate(contigs):
        y = len(contigs) - i - 1
        if colour_by in ("gc", "depth"):
            if colour_by == "gc":
                values = np.asarray(
                    [_gc(contig.sequence[s:e]) for s, e in _windows(contig.length, window)]
                )
                lo, hi = gc_range
            else:
                values = depths[contig.name]
                lo, hi = 0.0, float(np.percentile(values, 98)) or 1.0
            span = (hi - lo) or 1.0
            shades = cmap(np.clip((values - lo) / span, 0, 1))
            edges = np.asarray([s for s, _ in _windows(contig.length, window)], dtype=float)
            widths = np.asarray([e - s for s, e in _windows(contig.length, window)], dtype=float)
            ax.bar(
                edges,
                bar_h,
                width=widths,
                bottom=y - bar_h / 2,
                align="edge",
                color=shades,
                linewidth=0,
            )
        else:
            colour = (
                t["palette"][2]
                if (colour_by == "shape" and contig.is_circular)
                else t["palette"][0]
            )
            ax.bar(
                0,
                bar_h,
                width=contig.length,
                bottom=y - bar_h / 2,
                align="edge",
                color=colour,
                linewidth=0,
            )

        # An outline, heavier when the contig is one the caller asked to see.
        # A plain rectangle, because a rounded box rounds in data coordinates
        # and the x and y scales here differ by six orders of magnitude, which
        # turned every corner into a stray vertical line across the figure.
        picked = contig.name in highlight
        ax.add_patch(
            Rectangle(
                (0, y - bar_h / 2),
                contig.length,
                bar_h,
                linewidth=1.6 if picked else 0.7,
                edgecolor=t["palette"][1] if picked else t["grid"],
                facecolor="none",
                zorder=3,
            )
        )

        label = contig.name
        if contig.is_circular:
            label += "  ○"  # a small ring, for the ones that closed
        ax.text(
            -longest * 0.012,
            y,
            label,
            ha="right",
            va="center",
            fontsize=8,
            color=t["palette"][1] if picked else t["fg"],
            fontweight="bold" if picked else "normal",
        )
        ax.text(
            contig.length + longest * 0.012,
            y,
            bp_label(contig.length, longest),
            ha="left",
            va="center",
            fontsize=7.5,
            color=t["muted"],
        )

    ax.set_ylim(-0.8, len(contigs) - 0.2)
    ax.set_xlim(-longest * 0.22, longest * 1.16)
    ax.set_yticks([])
    for side in ("top", "right", "left"):
        ax.spines[side].set_visible(False)
    ax.set_xlabel("position")
    ticks = np.linspace(0, longest, 6)
    ax.set_xticks(ticks)
    ax.set_xticklabels([bp_label(v, longest) for v in ticks])
    ax.grid(axis="x", alpha=0.25, linewidth=0.6, color=t["grid"])
    ax.set_axisbelow(True)

    if colour_by in ("gc", "depth"):
        mappable = plt.cm.ScalarMappable(cmap=cmap)
        mappable.set_array([])
        bar = fig.colorbar(mappable, ax=ax, fraction=0.025, pad=0.02)
        bar.set_label("GC content" if colour_by == "gc" else "read depth", color=t["fg"])
        bar.ax.tick_params(colors=t["muted"])
        bar.set_ticks([0, 0.5, 1])
        if colour_by == "gc":
            lo, hi = gc_range
            bar.set_ticklabels([f"{lo:.0%}", f"{(lo + hi) / 2:.0%}", f"{hi:.0%}"])
        else:
            bar.set_ticklabels(["low", "", "high"])
        bar.outline.set_edgecolor(t["grid"])

    heading = title or f"{Path(fasta).stem}"
    subtitle = f"{n_all:,} contigs, {bp_label(total_all)} in total"
    if hidden:
        subtitle += f". The largest {len(contigs)} drawn, {hidden:,} smaller ones left out"
    ax.set_title(f"{heading}\n{subtitle}", fontsize=11, pad=14)

    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=dpi, bbox_inches="tight", facecolor=fig.get_facecolor())
    plt.close(fig)
    return path


def _depth_per_contig(
    contigs: list[Contig],
    reads: list[str] | None,
    window: int,
    read_type: str = "ont",
    threads: int = 4,
    max_reads: int = 20_000,
) -> dict[str, np.ndarray] | None:
    """Windowed depth for every contig at once, from one minimap2 run."""
    if not reads or not have("minimap2"):
        return None

    tmp = Path(tempfile.mkdtemp(prefix="nagasm_karyo_"))
    target, reads_path = tmp / "contigs.fasta", tmp / "reads.fastq"
    write_fasta(target, [Record(c.name, c.sequence) for c in contigs], wrap=0)

    factor = 1.0
    if len(reads) > max_reads:
        step = len(reads) / max_reads
        factor = step
        reads = [reads[int(i * step)] for i in range(max_reads)]
    write_fastq(reads_path, [Record(f"read{i}", r) for i, r in enumerate(reads)])

    proc = subprocess.run(
        ["minimap2", "-x", _preset_for(read_type), "-t", str(threads), str(target), str(reads_path)],
        stdout=subprocess.PIPE,
        stderr=subprocess.DEVNULL,
    )
    if proc.returncode != 0:
        return None

    out = {
        c.name: np.zeros(len(_windows(c.length, window)), dtype=float) for c in contigs
    }
    lengths = {c.name: c.length for c in contigs}
    for line in proc.stdout.decode().splitlines():
        fields = line.split("\t")
        if len(fields) < 11:
            continue
        name = fields[5]
        if name not in out:
            continue
        start, end = int(fields[7]), int(fields[8])
        first, last = start // window, min((end - 1) // window, len(out[name]) - 1)
        for w in range(first, last + 1):
            lo, hi = w * window, min((w + 1) * window, lengths[name])
            out[name][w] += max(0, min(end, hi) - max(start, lo))

    for name, values in out.items():
        spans = np.asarray([e - s for s, e in _windows(lengths[name], window)], dtype=float)
        out[name] = values * factor / np.maximum(spans, 1.0)
    return out


# ------------------------------------------------------------------ dotplot


def _offsets(fasta: Path) -> tuple[dict[str, int], int]:
    """Where each sequence starts when they are laid end to end."""
    offsets, running = {}, 0
    for record in read_sequences(fasta):
        offsets[record.name] = running
        running += len(record.seq)
    return offsets, running


def dotplot(
    query: Path,
    path: Path,
    reference: Path | None = None,
    min_block: int = 200,
    theme_name: str = "light",
    title: str = "",
    threads: int = 4,
    figsize: tuple[float, float] = (8.0, 8.0),
    dpi: int = 150,
    colour_by_identity: bool = True,
) -> Path | None:
    """A sequence against itself, or against a reference.

    Against itself the main diagonal is removed, so everything left is a
    repeat. That is the useful version, because repeats are what break
    assemblies and a dotplot shows you their length and their copy number at a
    glance.

    Returns None if minimap2 is missing or nothing aligned, because "nothing
    aligned" is a real answer.
    """
    if not have("minimap2"):
        return None

    self_mode = reference is None
    target = Path(query) if self_mode else Path(reference)
    cmd = ["minimap2", "-x", "asm20", "-t", str(threads)]
    if self_mode:
        # -D drops the self diagonal, -P keeps every chain rather than the best
        # one, which is exactly what you want when looking for repeats.
        cmd += ["-D", "-P"]
    cmd += [str(target), str(Path(query))]

    proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
    if proc.returncode != 0:
        return None

    # A FASTA with several sequences in it has to be laid end to end, with an
    # offset per sequence. Without that, position 500 of a 3 kb circle and
    # position 500 of a 48 kb one land on the same spot and the picture is a
    # lie. Separators are drawn between them.
    q_offsets, q_total = _offsets(Path(query))
    t_offsets, t_total = (q_offsets, q_total) if self_mode else _offsets(target)

    segments = []
    for line in proc.stdout.decode().splitlines():
        fields = line.split("\t")
        if len(fields) < 11:
            continue
        qname, tname = fields[0], fields[5]
        if qname not in q_offsets or tname not in t_offsets:
            continue
        qo, to = q_offsets[qname], t_offsets[tname]
        qs, qe = int(fields[2]) + qo, int(fields[3]) + qo
        strand = fields[4]
        ts, te = int(fields[7]) + to, int(fields[8]) + to
        block = int(fields[10])
        if block < min_block:
            continue
        identity = int(fields[9]) / block if block else 0.0
        segments.append((qs, qe, ts, te, strand, identity))

    if not segments:
        return None
    qlen, tlen = q_total, t_total

    fig, ax = plt.subplots(figsize=figsize)
    t = apply_theme(fig, ax, theme_name)
    cmap = gradient("identity")

    for qs, qe, ts, te, strand, identity in segments:
        xs = (qs, qe)
        ys = (ts, te) if strand == "+" else (te, ts)
        colour = cmap(np.clip((identity - 0.80) / 0.20, 0, 1)) if colour_by_identity else (
            t["palette"][0] if strand == "+" else t["palette"][1]
        )
        ax.plot(xs, ys, color=colour, linewidth=1.3, solid_capstyle="round")

    ax.set_xlim(0, qlen)
    ax.set_ylim(0, tlen)
    ax.set_xlabel(Path(query).stem)
    ax.set_ylabel(Path(query).stem if self_mode else Path(reference).stem)

    # Where one sequence ends and the next begins, with its name, so long as
    # there are few enough to label without a mess.
    def _mark(axis_set_ticks, axis_set_labels, offsets, total, vertical: bool):
        if len(offsets) <= 1:
            ticks = np.linspace(0, total, 6)
            axis_set_ticks(ticks)
            axis_set_labels([bp_label(v, total) for v in ticks])
            return
        starts = sorted(offsets.values())
        for edge in starts[1:]:
            if vertical:
                ax.axvline(edge, color=t["grid"], linewidth=0.8)
            else:
                ax.axhline(edge, color=t["grid"], linewidth=0.8)
        if len(offsets) <= 12:
            ordered = sorted(offsets.items(), key=lambda kv: kv[1])
            bounds = [o for _, o in ordered] + [total]
            centres = [(bounds[i] + bounds[i + 1]) / 2 for i in range(len(ordered))]
            axis_set_ticks(centres)
            axis_set_labels([n for n, _ in ordered], fontsize=7, rotation=0 if not vertical else 30)
        else:
            ticks = np.linspace(0, total, 6)
            axis_set_ticks(ticks)
            axis_set_labels([bp_label(v, total) for v in ticks])

    _mark(ax.set_xticks, ax.set_xticklabels, q_offsets, qlen, True)
    _mark(ax.set_yticks, ax.set_yticklabels, t_offsets, tlen, False)
    ax.grid(alpha=0.2, linewidth=0.6, color=t["grid"])
    ax.set_axisbelow(True)

    heading = title or ("self comparison, repeats only" if self_mode else "against reference")
    ax.set_title(f"{heading}\n{len(segments):,} alignment blocks over {min_block} bp", fontsize=11)

    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=dpi, bbox_inches="tight", facecolor=fig.get_facecolor())
    plt.close(fig)
    return path


# ----------------------------------------------------------------- coverage


def coverage_plot(
    fasta: Path,
    reads: list[str],
    path: Path,
    window: int | None = None,
    max_contigs: int = 8,
    theme_name: str = "light",
    title: str = "",
    read_type: str = "ont",
    threads: int = 4,
    dpi: int = 150,
    log_scale: bool = False,
) -> Path | None:
    """Read depth along each contig, drawn straight rather than round.

    A flat line is what a real replicon looks like. A step in the middle is
    usually two things joined that should not have been, and a spike is usually
    a collapsed repeat.
    """
    contigs = _load_contigs(Path(fasta))
    if not contigs or not reads:
        return None
    contigs.sort(key=lambda c: c.length, reverse=True)
    contigs = contigs[:max_contigs]
    longest = max(c.length for c in contigs)
    if window is None:
        window = max(50, longest // 600)

    depths = _depth_per_contig(contigs, reads, window, read_type=read_type, threads=threads)
    if depths is None:
        return None

    fig, axes = plt.subplots(
        len(contigs), 1, figsize=(10, max(2.0, 1.5 * len(contigs))), sharex=False
    )
    axes = np.atleast_1d(axes)
    t = apply_theme(fig, axes, theme_name)

    for ax, contig in zip(axes, contigs):
        values = depths[contig.name]
        x = np.asarray([s for s, _ in _windows(contig.length, window)], dtype=float)
        ax.fill_between(x, values, color=t["palette"][0], alpha=0.75, linewidth=0)
        median = float(np.median(values))
        ax.axhline(median, color=t["palette"][1], linewidth=1.0, linestyle="--")
        ax.text(
            contig.length,
            median,
            f"  median {median:.0f}x",
            va="center",
            fontsize=7.5,
            color=t["palette"][1],
        )
        if log_scale and values.max() > 0:
            ax.set_yscale("log")
        ax.set_xlim(0, contig.length)
        ax.set_ylabel("depth", fontsize=8)
        shape = " (circular)" if contig.is_circular else ""
        ax.set_title(
            f"{contig.name}{shape}, {bp_label(contig.length)}",
            fontsize=9,
            loc="left",
        )
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
        ticks = np.linspace(0, contig.length, 6)
        ax.set_xticks(ticks)
        ax.set_xticklabels([bp_label(v, contig.length) for v in ticks], fontsize=7.5)
        ax.grid(alpha=0.2, linewidth=0.6, color=t["grid"])
        ax.set_axisbelow(True)

    if title:
        fig.suptitle(title, color=t["fg"], fontsize=11)
    fig.tight_layout()
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=dpi, bbox_inches="tight", facecolor=fig.get_facecolor())
    plt.close(fig)
    return path


# ----------------------------------------------------------------------- nx


def nx_curve(
    fasta: Path,
    path: Path,
    compare: dict[str, Path] | None = None,
    theme_name: str = "light",
    title: str = "",
    figsize: tuple[float, float] = (8.0, 5.0),
    dpi: int = 150,
) -> Path:
    """The contiguity curve, of which N50 is a single point.

    Read it as: for each x along the bottom, the y is the length such that
    contigs at least that long make up x percent of the assembly. A curve that
    stays high for a long way is a good assembly. N50 is where it crosses 50,
    and quoting only that number hides whether the curve falls off a cliff just
    after it.
    """
    # Keys have to be unique or one assembly silently replaces another, and
    # two runs of this pipeline both produce a file called nuclear.fasta.
    sets: dict[str, Path] = {}
    for source in [Path(fasta)] + [Path(v) for v in (compare or {}).values()]:
        label = source.stem
        depth = 1
        # Walk back up the path until the label is unique. Two files both
        # called nuclear.fasta in a directory both called assembly need three
        # components before they differ.
        while label in sets and depth < 5:
            depth += 1
            label = "/".join(source.resolve().parts[-depth:]).removesuffix(source.suffix)
        n = 2
        while label in sets:
            label, n = f"{label} ({n})", n + 1
        sets[label] = source

    fig, ax = plt.subplots(figsize=figsize)
    t = apply_theme(fig, ax, theme_name)

    for i, (label, source) in enumerate(sets.items()):
        lengths = sorted((len(r.seq) for r in read_sequences(source)), reverse=True)
        if not lengths:
            continue
        total = sum(lengths)
        cumulative = np.cumsum(lengths) / total * 100
        colour = t["palette"][i % len(t["palette"])]
        # The curve has to start at zero percent. Plotting the cumulative
        # values alone begins the line at the first contig's own share, which
        # for a bacterial chromosome is 95 percent, leaving the left hand
        # ninety five percent of the plot blank.
        x = np.concatenate(([0.0], cumulative))
        y = np.concatenate(([lengths[0]], lengths))
        ax.step(x, y, where="post", color=colour, linewidth=1.8, label=label)

        n50 = next(l for c, l in zip(cumulative, lengths) if c >= 50)
        ax.plot([50], [n50], "o", color=colour, markersize=6, zorder=4)
        ax.annotate(
            f"N50 {bp_label(n50)}",
            (50, n50),
            textcoords="offset points",
            xytext=(8, 8),
            fontsize=8,
            color=colour,
        )

    ax.set_yscale("log")
    ax.set_xlim(0, 100)
    ax.set_xlabel("percent of the assembly, in contigs at least this long")
    ax.set_ylabel("contig length")
    ax.axvline(50, color=t["grid"], linewidth=0.8, linestyle=":")
    ax.grid(alpha=0.25, linewidth=0.6, color=t["grid"])
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    if len(sets) > 1:
        legend = ax.legend(frameon=False, fontsize=9)
        for text in legend.get_texts():
            text.set_color(t["fg"])
    ax.set_title(title or "contiguity", fontsize=11, color=t["fg"])

    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=dpi, bbox_inches="tight", facecolor=fig.get_facecolor())
    plt.close(fig)
    return path
