"""Drawing a circular genome as a circle.

This is an optional extra. Nothing in the assembly depends on it, and the
pipeline runs fine if you never call it. It exists because a closed circle
drawn as a straight line hides the one thing that makes it interesting, which
is that the two ends meet.

There are two modes, and the difference between them is worth understanding
before you pick one.

**Reference free.** Everything drawn comes from the element itself and, if you
hand them over, the reads. GC content, GC skew and read depth around the
circle. This is the honest mode, because a new isolate has no reference to
compare against and this still works. GC skew is the interesting track for a
real replicon, since a bacterial chromosome usually flips sign at the origin
and the terminus, giving two clean halves. A circle that is actually a
collapsed repeat tends to show nothing of the sort.

**Reference based.** Adds rings built by aligning the element to a reference
you supply. One ring for where the alignment lands, shaded by identity, and
one for the gaps. This answers a different question: not is my circle real,
but how does it differ from the one already published. It needs minimap2 and
it needs you to already have the answer, which is exactly the dependency the
rest of this pipeline is built to avoid. Useful for a write up, not for
discovery.

A note on what the picture can and cannot tell you. The rings are drawn from
the same alignments the support check uses, so a plot that looks beautiful is
not independent evidence that the element is real. Read the containment number
next to it. The plot is for showing someone else what you found, and for
spotting the odd case where a number looks fine and the structure does not.
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

from .seqio import Record, read_sequences, write_fasta, write_fastq
from .polish import have, _preset_for
from .plots import PALETTE as _COLOURS, apply_theme, gradient, bp_label, safe_name

_IDENTITY_CMAP = gradient("identity")

# The rings that can be drawn, in the order they are drawn from the outside in.
# Pass a subset to plot_circle to leave some out.
ALL_RINGS = ("gc", "skew", "cumulative", "depth", "identity", "missing")


@dataclass
class Ring:
    """One track to draw around the circle.

    `values` is one number per window, already in the order the windows run
    around the circle. `vmin` and `vmax` fix the colour or height scale so that
    two plots of different elements can be put side by side.
    """

    label: str
    values: np.ndarray
    kind: str = "deviation"  # deviation, curve, heat or blocks
    colour: str = _COLOURS[0]
    vmin: float | None = None
    vmax: float | None = None
    # Where the baseline sits for a deviation ring. None means use the mean,
    # which is what you want for GC content, since the question is always
    # "unusual for this element" rather than "unusual in absolute terms".
    centre: float | None = None
    note: str = ""


def _windows(length: int, window: int) -> list[tuple[int, int]]:
    if window <= 0:
        raise ValueError("window must be positive")
    return [(s, min(s + window, length)) for s in range(0, length, window)]


def gc_content(sequence: str, window: int) -> np.ndarray:
    """GC fraction in each window, as a fraction not a percentage."""
    out = []
    for start, end in _windows(len(sequence), window):
        chunk = sequence[start:end].upper()
        gc = chunk.count("G") + chunk.count("C")
        at = chunk.count("A") + chunk.count("T")
        out.append(gc / (gc + at) if gc + at else 0.0)
    return np.asarray(out, dtype=float)


def gc_skew(sequence: str, window: int) -> np.ndarray:
    """(G - C) / (G + C) in each window.

    On a real bacterial chromosome this changes sign at the replication origin
    and again at the terminus, because the leading strand accumulates G. Two
    clean halves is a strong hint that you are looking at a replicon that is
    actually replicated. It is not proof, and plenty of real elements show
    nothing, so a flat skew is not evidence against anything.
    """
    out = []
    for start, end in _windows(len(sequence), window):
        chunk = sequence[start:end].upper()
        g, c = chunk.count("G"), chunk.count("C")
        out.append((g - c) / (g + c) if g + c else 0.0)
    return np.asarray(out, dtype=float)


def cumulative_gc_skew(sequence: str, window: int) -> np.ndarray:
    """Running total of the skew, which is where the origin shows up clearly.

    The raw skew is noisy. Its running total is not, and its minimum and
    maximum sit near the origin and the terminus.
    """
    return np.cumsum(gc_skew(sequence, window))


def read_depth(
    sequence: str,
    reads: list[str],
    window: int,
    is_circular: bool = True,
    read_type: str = "ont",
    threads: int = 4,
    max_reads: int = 8_000,
    workdir: Path | None = None,
) -> np.ndarray | None:
    """Depth in each window, from mapping the reads back.

    Returns None if minimap2 is missing, so the plot drops the ring rather than
    failing. The element is doubled for the same reason the support check
    doubles it, which is that a read crossing the join needs somewhere to land.
    Alignments are then folded back onto the real coordinates.
    """
    if not reads or not have("minimap2"):
        return None

    true_length = len(sequence)
    target = sequence + sequence if is_circular else sequence

    factor = 1.0
    if len(reads) > max_reads:
        step = len(reads) / max_reads
        factor = step
        reads = [reads[int(i * step)] for i in range(max_reads)]

    tmp = Path(workdir) if workdir else Path(tempfile.mkdtemp(prefix="nagasm_circle_"))
    tmp.mkdir(parents=True, exist_ok=True)
    target_path, reads_path = tmp / "target.fasta", tmp / "reads.fastq"
    write_fasta(target_path, [Record("element", target)], wrap=0)
    write_fastq(reads_path, [Record(f"read{i}", r) for i, r in enumerate(reads)])

    proc = subprocess.run(
        [
            "minimap2",
            "-x",
            _preset_for(read_type),
            "-t",
            str(threads),
            str(target_path),
            str(reads_path),
        ],
        stdout=subprocess.PIPE,
        stderr=subprocess.DEVNULL,
    )
    if proc.returncode != 0:
        return None

    bounds = _windows(true_length, window)
    covered = np.zeros(len(bounds), dtype=float)
    for line in proc.stdout.decode().splitlines():
        fields = line.split("\t")
        if len(fields) < 11:
            continue
        start, end = int(fields[7]), int(fields[8])
        for pos in range(start, end):
            # Fold the doubled coordinate back onto the real circle.
            covered[min((pos % true_length) // window, len(bounds) - 1)] += 1.0

    spans = np.asarray([e - s for s, e in bounds], dtype=float)
    return covered * factor / np.maximum(spans, 1.0)


def reference_rings(
    sequence: str,
    reference: Path,
    window: int,
    threads: int = 4,
    workdir: Path | None = None,
) -> list[Ring]:
    """Rings showing how the element lines up against a reference.

    Two rings come back: identity where the element aligns, and a flag for the
    windows where it does not align at all. Both are empty if minimap2 is
    missing or nothing aligns, and an empty list is a valid answer rather than
    an error, because "this does not match your reference" is information.
    """
    if not have("minimap2"):
        return []

    tmp = Path(workdir) if workdir else Path(tempfile.mkdtemp(prefix="nagasm_circleref_"))
    tmp.mkdir(parents=True, exist_ok=True)
    query = tmp / "element.fasta"
    write_fasta(query, [Record("element", sequence)], wrap=0)

    proc = subprocess.run(
        ["minimap2", "-x", "asm20", "-t", str(threads), str(reference), str(query)],
        stdout=subprocess.PIPE,
        stderr=subprocess.DEVNULL,
    )
    if proc.returncode != 0:
        return []

    bounds = _windows(len(sequence), window)
    best = np.zeros(len(bounds), dtype=float)
    for line in proc.stdout.decode().splitlines():
        fields = line.split("\t")
        if len(fields) < 11:
            continue
        qstart, qend = int(fields[2]), int(fields[3])
        matches, block = int(fields[9]), int(fields[10])
        identity = matches / block if block else 0.0
        for i, (s, e) in enumerate(bounds):
            if s < qend and e > qstart:
                best[i] = max(best[i], identity)

    if not best.any():
        return []

    missing = (best == 0).astype(float)
    hit = best[best > 0]
    note = f"{hit.min():.1%} to {hit.max():.1%}"
    rings = [
        Ring("identity to reference", best, kind="heat", vmin=0.80, vmax=1.0, note=note),
    ]
    if missing.any():
        share = missing.mean()
        rings.append(
            Ring(
                "not in reference",
                missing,
                kind="blocks",
                colour=_COLOURS[1],
                note=f"{share:.0%} of the element",
            )
        )
    return rings


def _draw_ring(ax, ring: Ring, inner: float, outer: float) -> None:
    n = len(ring.values)
    if n == 0:
        return
    edges = np.linspace(0, 2 * math.pi, n + 1)
    width = edges[1] - edges[0]
    mid = edges[:-1] + width / 2

    if ring.kind == "heat":
        lo = ring.vmin if ring.vmin is not None else float(np.min(ring.values))
        hi = ring.vmax if ring.vmax is not None else float(np.max(ring.values))
        span = hi - lo or 1.0
        shades = _IDENTITY_CMAP(np.clip((ring.values - lo) / span, 0, 1))
        ax.bar(
            mid, outer - inner, width=width, bottom=inner, color=shades, linewidth=0
        )
        return

    if ring.kind == "blocks":
        on = ring.values > 0
        if on.any():
            ax.bar(
                mid[on],
                outer - inner,
                width=width,
                bottom=inner,
                color=ring.colour,
                linewidth=0,
            )
        return

    if ring.kind == "curve":
        # A continuous line, closed back on itself. This is the right shape for
        # anything cumulative, where what matters is where the peak and the
        # trough sit rather than the value in any one window.
        lo = float(np.min(ring.values))
        hi = float(np.max(ring.values))
        span = hi - lo or 1.0
        scaled = (ring.values - lo) / span
        radius = inner + scaled * (outer - inner)
        closed_angle = np.append(mid, mid[0] + 2 * math.pi)
        closed_radius = np.append(radius, radius[0])
        ax.plot(closed_angle, closed_radius, color=ring.colour, linewidth=1.2, zorder=2)
        # Mark the extremes. On a bacterial chromosome these sit near the
        # origin and the terminus of replication.
        for idx, marker in ((int(np.argmax(ring.values)), "^"), (int(np.argmin(ring.values)), "v")):
            ax.plot(
                mid[idx], radius[idx], marker, color=ring.colour, markersize=5, zorder=3
            )
        return

    # A deviation ring: bars growing out from a baseline, so that above and
    # below the baseline are visibly different things rather than just taller
    # and shorter bars.
    values = ring.values
    centre = float(np.mean(values)) if ring.centre is None else ring.centre
    reach = ring.vmax if ring.vmax is not None else float(np.max(np.abs(values - centre)))
    reach = reach or 1.0
    scaled = np.clip((values - centre) / reach, -1, 1)
    mid_radius = (inner + outer) / 2
    half = (outer - inner) / 2
    ax.plot(
        np.linspace(0, 2 * math.pi, 200),
        np.full(200, mid_radius),
        color="#BBBBBB",
        linewidth=0.5,
        zorder=1,
    )
    ax.bar(
        mid,
        scaled * half,
        width=width,
        bottom=mid_radius,
        color=ring.colour,
        linewidth=0,
        zorder=2,
    )


def plot_circle(
    sequence: str,
    path: Path,
    name: str = "element",
    reads: list[str] | None = None,
    reference: Path | None = None,
    window: int | None = None,
    is_circular: bool = True,
    read_type: str = "ont",
    threads: int = 4,
    note: str = "",
    rings: tuple[str, ...] | list[str] | None = None,
    theme_name: str = "light",
    title: str = "",
    legend: bool = True,
    start: int = 0,
    figsize: tuple[float, float] = (8.5, 8.5),
    dpi: int = 150,
) -> Path:
    """Draw one circular element and save it.

    Rings go from the outside in: GC content, GC skew, cumulative skew, then
    read depth if reads were given, then the reference rings if a reference was
    given. A ring whose data could not be built is left out rather than drawn
    empty.
    """
    length = len(sequence)
    if length < 100:
        raise ValueError("too short to draw a useful circle")
    if start:
        # Rotate the circle so a feature you care about sits at the top. Only
        # meaningful for something actually circular.
        start %= length
        sequence = sequence[start:] + sequence[:start]
    wanted = set(rings) if rings else set(ALL_RINGS)
    if window is None:
        # Aim for a few hundred windows, which looks right at any size from a
        # 3 kb circle to a 5 Mb chromosome.
        window = max(10, length // 360)

    palette = _theme_palette(theme_name)
    gc = gc_content(sequence, window)
    rings: list[Ring] = []
    if "gc" in wanted:
        rings.append(
            Ring("GC content", gc, colour=palette[0], note=f"mean {gc.mean():.1%}")
        )
    if "skew" in wanted:
        rings.append(
            Ring(
                "GC skew",
                gc_skew(sequence, window),
                colour=palette[2],
                centre=0.0,
                note="(G-C)/(G+C), around zero",
            )
        )
    if "cumulative" in wanted:
        rings.append(
            Ring(
                "cumulative GC skew",
                cumulative_gc_skew(sequence, window),
                kind="curve",
                colour=palette[5],
                note="peak and trough marked",
            )
        )

    if reads and "depth" in wanted:
        depth = read_depth(
            sequence,
            reads,
            window,
            is_circular=is_circular,
            read_type=read_type,
            threads=threads,
        )
        if depth is not None:
            rings.append(
                Ring(
                    "read depth",
                    depth,
                    colour=palette[4],
                    note=f"mean {depth.mean():.0f}x, around the mean",
                )
            )

    if reference is not None:
        for ring in reference_rings(sequence, Path(reference), window, threads=threads):
            key = "identity" if ring.kind == "heat" else "missing"
            if key in wanted:
                rings.append(ring)

    if not rings:
        raise ValueError(f"no rings left to draw. Pick from {', '.join(ALL_RINGS)}")

    fig = plt.figure(figsize=figsize)
    ax = fig.add_subplot(111, projection="polar")
    t = apply_theme(fig, ax, theme_name)
    ax.set_theta_zero_location("N")
    ax.set_theta_direction(-1)

    # Rings laid out from the outside in, with a gap between each.
    top, floor, gap = 1.0, 0.42, 0.015
    band = (top - floor) / len(rings)
    labels = []
    for i, ring in enumerate(rings):
        outer = top - i * band
        inner = outer - band + gap
        _draw_ring(ax, ring, inner, outer)
        labels.append((ring, (inner + outer) / 2))

    ax.set_ylim(0, 1.02)
    ax.set_yticks([])
    ax.grid(False)
    ax.spines["polar"].set_visible(False)

    ticks = np.linspace(0, 2 * math.pi, 9)[:-1]
    ax.set_xticks(ticks)
    ax.set_xticklabels(
        [bp_label(int(a / (2 * math.pi) * length), length) for a in ticks],
        fontsize=8,
        color=t["muted"],
    )

    # Labels go in a legend outside the circle. Written on the rings they sat
    # on top of the data they were describing.
    handles = []
    for ring, _ in labels:
        colour = ring.colour
        if ring.kind == "heat":
            colour = _IDENTITY_CMAP(0.75)
        text = ring.label if not ring.note else f"{ring.label}  ({ring.note})"
        handles.append(plt.Line2D([], [], color=colour, linewidth=5, label=text))
    if legend:
        drawn = ax.legend(
            handles=handles,
            loc="upper left",
            bbox_to_anchor=(-0.12, 1.06),
            frameon=False,
            fontsize=8,
            labelspacing=0.9,
            handlelength=1.4,
        )
        for text_item in drawn.get_texts():
            text_item.set_color(t["fg"])

    shape = "circular" if is_circular else "linear"
    heading = title or name
    caption = f"{heading}\n{length:,} bp, {shape}, {window:,} bp windows"
    if start:
        caption += f", rotated to start at {start:,}"
    # A rejected element is still drawn, because it is still reported, but the
    # whole title turns the warning colour so nobody mistakes it for one that
    # passed. Colouring only the warning line meant placing a second text box
    # by hand, and it landed on top of the title.
    failed = note.upper().startswith("NOT SUPPORTED")
    if note:
        caption += f"\n{note}"
    ax.set_title(
        caption,
        fontsize=11,
        pad=24,
        color=_COLOURS[1] if failed else t["fg"],
        fontweight="bold" if failed else "normal",
    )

    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=dpi, bbox_inches="tight", facecolor=fig.get_facecolor())
    plt.close(fig)
    return path


def _theme_palette(name: str) -> list[str]:
    from .plots import theme as _theme

    return _theme(name)["palette"]


def plot_all(
    elements: Path,
    outdir: Path,
    reads: list[str] | None = None,
    reference: Path | None = None,
    read_type: str = "ont",
    threads: int = 4,
    window: int | None = None,
    rings: tuple[str, ...] | list[str] | None = None,
    theme_name: str = "light",
    legend: bool = True,
    image_format: str = "png",
    dpi: int = 150,
    only_circular: bool = False,
    min_length: int = 100,
    notes: dict[str, str] | None = None,
) -> list[Path]:
    """Draw every element in a FASTA file, one picture each.

    `notes` carries a line to print under the title, keyed by element name.
    The pipeline uses it to say on the face of the picture when an element
    failed the read support check, because a rejected artefact still gets
    drawn and a beautiful circle with nothing written on it is exactly how
    someone ends up believing in one.
    """
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    written = []
    for record in read_sequences(elements):
        if len(record.seq) < max(100, min_length):
            continue
        circular = "circular" in record.name.lower()
        if only_circular and not circular:
            continue
        written.append(
            plot_circle(
                record.seq,
                outdir / f"circle_{safe_name(record.name)}.{image_format}",
                name=record.name,
                reads=reads,
                reference=reference,
                window=window,
                is_circular=circular,
                read_type=read_type,
                threads=threads,
                rings=rings,
                theme_name=theme_name,
                legend=legend,
                dpi=dpi,
                note=(notes or {}).get(record.name, ""),
            )
        )
    return written



