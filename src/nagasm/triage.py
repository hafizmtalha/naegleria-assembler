"""Sorting reads into compartments before anything is assembled.

This is the part of the package that exists because of Naegleria specifically.

A Naegleria cell carries its ribosomal RNA genes on a small circle of DNA rather
than on its chromosomes, and it holds something like four thousand copies of that
circle. The mitochondrial genome is a second circle, also in many copies. Neither
is big. Together they are a fraction of a percent of the length of the nuclear
genome. But because of their copy number they can account for a large share of
the reads, and in some libraries more reads than the chromosomes do.

Hand that pile to an ordinary assembler and it struggles, because almost
everything an assembler does assumes coverage is roughly even. The team behind
the 2019 nanopore assembly hit this and worked around it by mapping reads to a
known rDNA circle from a related species and removing whatever matched. That
works when you already have a reference. It does not help you with a new isolate.

So we find the compartments in the reads themselves.

## How

Count k-mers across the whole read set and give each read the median count of its
own k-mers. Reads from a four thousand copy circle are built from k-mers that are
seen far more often than k-mers from a single copy chromosome, so that one number
already separates them.

To decide where the boundaries go we look at the k-mer count histogram, which is
how many distinct k-mers are seen exactly once, twice, and so on. Each
compartment makes its own hump. Crucially, the area under a hump tells you how
much distinct sequence sits at that depth. The nuclear hump covers tens of
megabases. The rDNA hump covers around fifteen kilobases at a thousand times the
depth.

That area is what lets us tell the compartments apart without knowing anything
about the organism. The nuclear compartment is simply the one that spans the most
distinct sequence. It is not the tallest hump and it is not the one with the most
reads, and assuming either of those would get Naegleria wrong.
"""

from __future__ import annotations

from dataclasses import dataclass, field, asdict
from typing import Sequence

import numpy as np

from . import nagcore

# Sizes of the small circles that have been published for Naegleria. These are
# only used to add a note next to a compartment we found on our own. Nothing is
# ever classified because of them, and a compartment that matches nothing still
# gets reported.
_SIZE_HINTS: tuple[tuple[int, int, str], ...] = (
    (12_000, 20_000, "size is in the range published for the Naegleria rDNA circle"),
    (40_000, 60_000, "size is in the range published for the Naegleria mitochondrial genome"),
)


@dataclass
class Compartment:
    """One hump in the k-mer depth histogram, and what we make of it."""

    label: str
    peak_depth: float
    low_depth: float
    high_depth: float
    span_bp: int
    n_reads: int = 0
    total_bases: int = 0
    median_gc: float = 0.0
    median_read_length: float = 0.0
    hint: str = ""

    @property
    def depth_relative_to_nuclear(self) -> float:
        """Filled in by the triage step once the nuclear peak is known."""
        return getattr(self, "_rel", float("nan"))

    def as_dict(self) -> dict:
        d = asdict(self)
        d["depth_relative_to_nuclear"] = round(self.depth_relative_to_nuclear, 2)
        return d


@dataclass
class TriageResult:
    """Everything the triage step worked out, including its inputs."""

    compartments: list[Compartment]
    assignment: np.ndarray  # one entry per read: index into compartments, or -1
    read_depth: np.ndarray
    read_gc: np.ndarray
    read_length: np.ndarray
    read_support: np.ndarray  # k-mers per read that survived the depth floor
    histogram: np.ndarray
    error_threshold: int
    k: int
    scale: int
    hpc: bool = False
    nuclear_index: int = -1
    warnings: list[str] = field(default_factory=list)

    @property
    def nuclear(self) -> Compartment | None:
        if self.nuclear_index < 0:
            return None
        return self.compartments[self.nuclear_index]

    @property
    def high_copy(self) -> list[Compartment]:
        """Compartments sitting above the nuclear one, deepest first."""
        return [
            c
            for i, c in enumerate(self.compartments)
            if i != self.nuclear_index and c.peak_depth > self.nuclear_peak
        ]

    @property
    def nuclear_peak(self) -> float:
        nuc = self.nuclear
        return nuc.peak_depth if nuc else float("nan")

    def read_indices(self, label: str) -> np.ndarray:
        """Which reads landed in the compartment with this label."""
        for i, c in enumerate(self.compartments):
            if c.label == label:
                return np.flatnonzero(self.assignment == i)
        raise KeyError(f"no compartment called {label!r}")


def find_error_threshold(histogram: np.ndarray, max_search: int = 100) -> int:
    """Where the error k-mers stop and the real ones begin.

    A sequencing error creates k-mers that exist nowhere else, so the histogram
    starts with an enormous spike at a count of one that falls away steeply. The
    real sequence makes humps further right. The dip between the two is the
    sensible place to cut.

    We walk right from a count of one until the histogram stops falling. If it
    never stops, which happens when coverage is too low for a real hump to form,
    we give up and return 1 rather than inventing a cut.
    """
    hist = np.asarray(histogram, dtype=float)
    if hist.size < 3:
        return 1
    limit = min(max_search, hist.size - 1)
    for d in range(1, limit):
        if hist[d + 1] > hist[d]:
            return d
    return 1


def _log_profile(
    histogram: np.ndarray, low: int, n_bins: int = 160
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Rebin the histogram onto a log spaced depth axis.

    Depth spans several orders of magnitude here, from a single copy chromosome
    to a four thousand copy circle. On a linear axis the nuclear hump is squashed
    into the first pixel. On a log axis every compartment gets similar room, which
    is what makes one peak finder work for all of them.
    """
    hist = np.asarray(histogram, dtype=float)
    top = hist.size - 1
    if top <= low:
        return np.array([]), np.array([]), np.array([])
    edges = np.geomspace(max(low, 1), top, n_bins + 1)
    depths = np.arange(hist.size)
    keep = depths >= low
    idx = np.clip(np.digitize(depths[keep], edges) - 1, 0, n_bins - 1)
    weight = np.bincount(idx, weights=hist[keep], minlength=n_bins)
    centers = np.sqrt(edges[:-1] * edges[1:])
    return centers, weight, edges


def _find_peaks(profile: np.ndarray, smoothing: float, min_prominence: float) -> list[int]:
    """Peak positions in the smoothed log profile.

    The peaks are found on the log of the profile, not the profile itself, and
    that detail matters more here than anywhere else in the package.

    The compartments differ in span by orders of magnitude. In a real Naegleria
    library the nuclear compartment covers tens of megabases while the rDNA
    circle covers about fifteen kilobases, so the rDNA hump is well under a
    percent of the height of the nuclear one. Any rule of the form "ignore peaks
    shorter than some fraction of the tallest" deletes exactly the compartment we
    built this tool to find.

    On a log axis the question becomes "does this stand out from its
    surroundings", which is the right question and has nothing to do with how
    tall the biggest peak is. `min_prominence` is in log10 units, so 0.15 means a
    peak has to rise about one and a half fold above the dips on either side.

    scipy is used when it is there. The fallback is a plain local maximum test,
    so the package still works without it.
    """
    if profile.size == 0:
        return []
    log_profile = np.log10(1.0 + np.asarray(profile, dtype=float))
    smooth = _smooth(log_profile, smoothing)
    try:
        from scipy.signal import find_peaks

        peaks, _ = find_peaks(smooth, prominence=min_prominence)
        return list(peaks)
    except ImportError:
        out = []
        for i in range(1, smooth.size - 1):
            if smooth[i] >= smooth[i - 1] and smooth[i] > smooth[i + 1]:
                left = smooth[:i].min() if i else smooth[i]
                right = smooth[i + 1 :].min() if i + 1 < smooth.size else smooth[i]
                if smooth[i] - max(left, right) >= min_prominence:
                    out.append(i)
        return out


def _smooth(values: np.ndarray, sigma: float) -> np.ndarray:
    """Gaussian blur along one axis, with edges held rather than wrapped."""
    if sigma <= 0:
        return values.astype(float)
    try:
        from scipy.ndimage import gaussian_filter1d

        return gaussian_filter1d(values.astype(float), sigma, mode="nearest")
    except ImportError:
        radius = max(1, int(round(3 * sigma)))
        x = np.arange(-radius, radius + 1)
        kernel = np.exp(-0.5 * (x / sigma) ** 2)
        kernel /= kernel.sum()
        padded = np.pad(values.astype(float), radius, mode="edge")
        return np.convolve(padded, kernel, mode="valid")


def find_compartments(
    histogram: np.ndarray,
    scale: int,
    error_threshold: int | None = None,
    smoothing: float = 2.0,
    min_prominence: float = 0.15,
    min_span_bp: int = 1_000,
) -> tuple[list[Compartment], int]:
    """Pull the humps out of a k-mer depth histogram.

    Gives back the compartments, deepest last, and the index of the one taken to
    be nuclear.

    The nuclear compartment is the one covering the most distinct sequence, not
    the one with the tallest peak. For most organisms those are the same hump. For
    Naegleria they are not, and that difference is the whole reason this function
    exists.
    """
    hist = np.asarray(histogram, dtype=float)
    low = error_threshold if error_threshold is not None else find_error_threshold(hist)
    centers, profile, edges = _log_profile(hist, low)
    if centers.size == 0:
        return [], -1

    peak_bins = _find_peaks(profile, smoothing, min_prominence)
    if not peak_bins:
        # No structure to speak of. Treat everything above the error cut as one
        # compartment rather than reporting nothing.
        span = int(hist[low:].sum() * scale)
        weights = hist[low:]
        depths = np.arange(low, hist.size)
        peak = float(depths[np.argmax(weights)]) if weights.size else float(low)
        c = Compartment("nuclear", peak, float(low), float(hist.size - 1), span)
        return [c], 0

    # Boundaries sit at the lowest point between neighbouring peaks.
    smooth = _smooth(np.log10(1.0 + profile), smoothing)
    bounds: list[tuple[int, int]] = []
    for i, p in enumerate(peak_bins):
        left = 0 if i == 0 else peak_bins[i - 1] + int(np.argmin(smooth[peak_bins[i - 1] : p]))
        right = (
            profile.size
            if i == len(peak_bins) - 1
            else p + int(np.argmin(smooth[p : peak_bins[i + 1]]))
        )
        bounds.append((left, right))

    comps: list[Compartment] = []
    for (left, right), p in zip(bounds, peak_bins):
        d_lo = float(edges[left])
        d_hi = float(edges[min(right, len(edges) - 1)])
        i_lo = max(low, int(np.floor(d_lo)))
        i_hi = min(hist.size, int(np.ceil(d_hi)) + 1)
        span = int(hist[i_lo:i_hi].sum() * scale)
        comps.append(
            Compartment(
                label="",
                peak_depth=float(centers[p]),
                low_depth=d_lo,
                high_depth=d_hi,
                span_bp=span,
            )
        )

    # A hump covering only a few hundred bases is noise, not a compartment.
    # Anything at or above min_span_bp is kept, however short it is next to the
    # nuclear compartment.
    comps = [c for c in comps if c.span_bp >= min_span_bp]
    if not comps:
        return [], -1

    nuclear_index = int(np.argmax([c.span_bp for c in comps]))
    _name_compartments(comps, nuclear_index)
    return comps, nuclear_index


def _name_compartments(comps: list[Compartment], nuclear_index: int) -> None:
    """Give each compartment a label and, where it fits, a size note."""
    nuclear_peak = comps[nuclear_index].peak_depth
    above = 0
    below = 0
    for i, c in enumerate(comps):
        c._rel = c.peak_depth / nuclear_peak if nuclear_peak else float("nan")
        if i == nuclear_index:
            c.label = "nuclear"
            continue
        if c.peak_depth > nuclear_peak:
            above += 1
            c.label = f"high_copy_{above}"
            for lo, hi, note in _SIZE_HINTS:
                if lo <= c.span_bp <= hi:
                    c.hint = note
                    break
        else:
            below += 1
            c.label = f"low_copy_{below}"
            c.hint = (
                "below the nuclear depth, which can mean a partial duplication, "
                "a second haplotype or something from another organism"
            )


def triage(
    reads: Sequence[str],
    k: int = 21,
    scale: int = 100,
    smoothing: float = 2.0,
    min_prominence: float = 0.15,
    max_depth: int = 20_000,
    batch_size: int = 200_000,
    hpc: bool | None = None,
    min_support: int = 3,
) -> TriageResult:
    """Sort reads into compartments using nothing but the reads.

    `scale` keeps roughly one k-mer in that many. It is the knob that decides how
    much memory this takes. The default keeps a hundredth of them, which is still
    tens of k-mers per long read and far more than a median needs.

    `max_depth` only sets where the histogram is truncated. Anything deeper is
    piled into the last slot, so a very deep compartment is still seen.

    `hpc` collapses runs of the same base before counting. Long read platforms
    lose count inside homopolymers, so the same stretch of DNA read twice gives
    two different k-mers and neither looks shared. Collapsing the runs first makes
    that disagreement go away. It is picked automatically when it is not set:
    on if the reads look long and noisy, off otherwise.

    `min_support` is how many k-mers a read needs above the error cut before its
    depth is believed. Reads below it are left unassigned rather than guessed at.
    """
    reads = list(reads)
    if not reads:
        raise ValueError("no reads given")

    if hpc is None:
        hpc = _looks_like_noisy_long_reads(reads)
    counted = [nagcore.homopolymer_compress(r)[0] for r in reads] if hpc else reads

    index = nagcore.DepthIndex(k=k, scale=scale)
    for start in range(0, len(counted), batch_size):
        index.add_reads(counted[start : start + batch_size])

    hist = np.asarray(index.depth_histogram(max_depth), dtype=np.int64)
    error_threshold = find_error_threshold(hist)
    comps, nuclear_index = find_compartments(
        hist, scale, error_threshold, smoothing, min_prominence
    )

    # Profile with the error cut as a floor. Without it, a noisy read's median is
    # the depth of its own sequencing errors rather than the depth of the
    # compartment it came from. See DepthIndex.profile_into for the detail.
    profiles = index.profile_reads(counted, min_count=max(error_threshold, 1))
    depth = np.array([p[0] for p in profiles], dtype=float)
    gc = np.array([p[2] for p in profiles], dtype=float)
    length = np.array([len(r) for r in reads], dtype=np.int64)
    support = np.array([p[5] for p in profiles], dtype=np.int64)

    assignment = _assign_reads(depth, support, comps, min_support)
    _fill_in_read_stats(comps, assignment, depth, gc, length)

    warnings = _sanity_checks(comps, nuclear_index, assignment, support, k, scale, hpc)

    return TriageResult(
        compartments=comps,
        assignment=assignment,
        read_depth=depth,
        read_gc=gc,
        read_length=length,
        read_support=support,
        histogram=hist,
        error_threshold=error_threshold,
        k=k,
        scale=scale,
        hpc=hpc,
        nuclear_index=nuclear_index,
        warnings=warnings,
    )


def _looks_like_noisy_long_reads(reads: Sequence[str], sample: int = 2000) -> bool:
    """Guess whether these are long reads that would benefit from collapsing.

    Short accurate reads do not need it and lose a little real signal from it, so
    the test is simply whether the reads are long. Anything with a typical length
    over a kilobase came off a long read platform.
    """
    lengths = [len(r) for r in reads[:sample]]
    if not lengths:
        return False
    return float(np.median(lengths)) >= 1000


def _assign_reads(
    depth: np.ndarray,
    support: np.ndarray,
    comps: list[Compartment],
    min_support: int,
) -> np.ndarray:
    """Send each read to the compartment whose peak is nearest in log depth.

    Distance is measured in log space because the question that matters is "how
    many times deeper", not "how many more copies". A read at depth 40 belongs with
    a peak at 30 rather than one at 900, and plain subtraction would get that
    wrong once the peaks are far apart.
    """
    assignment = np.full(depth.size, -1, dtype=np.int64)
    if not comps:
        return assignment

    peaks = np.array([c.peak_depth for c in comps], dtype=float)
    log_peaks = np.log(np.maximum(peaks, 1e-9))

    usable = (support >= max(min_support, 1)) & (depth > 0)
    if not usable.any():
        return assignment

    log_depth = np.log(np.maximum(depth[usable], 1e-9))
    nearest = np.argmin(np.abs(log_depth[:, None] - log_peaks[None, :]), axis=1)
    assignment[usable] = nearest
    return assignment


def _fill_in_read_stats(
    comps: list[Compartment],
    assignment: np.ndarray,
    depth: np.ndarray,
    gc: np.ndarray,
    length: np.ndarray,
) -> None:
    for i, c in enumerate(comps):
        sel = assignment == i
        c.n_reads = int(sel.sum())
        c.total_bases = int(length[sel].sum())
        c.median_gc = float(np.median(gc[sel])) if c.n_reads else 0.0
        c.median_read_length = float(np.median(length[sel])) if c.n_reads else 0.0


def _sanity_checks(
    comps: list[Compartment],
    nuclear_index: int,
    assignment: np.ndarray,
    support: np.ndarray,
    k: int,
    scale: int,
    hpc: bool,
) -> list[str]:
    """Things worth telling the user about rather than silently carrying on."""
    notes: list[str] = []

    unassigned = int((assignment < 0).sum())
    if unassigned:
        share = 100.0 * unassigned / assignment.size
        notes.append(
            f"{unassigned} reads ({share:.1f}%) had too few k-mers above the error "
            f"cut to place, so they are in no compartment."
        )
        if share > 20:
            advice = (
                f"That is a lot. It usually means the reads are noisier than k={k} "
                f"can cope with, because an error anywhere in a k-mer destroys it. "
                f"Try a smaller k, or a smaller scale than {scale} so each read "
                f"contributes more k-mers."
            )
            if not hpc:
                advice += " Turning on hpc would also help if these are long reads."
            notes.append(advice)

    if nuclear_index >= 0:
        nuc = comps[nuclear_index]
        if nuc.span_bp < 1_000_000:
            notes.append(
                f"The compartment taken to be nuclear only spans about "
                f"{nuc.span_bp / 1000:.0f} kb. That is small for a nuclear genome, so "
                f"check the depth plot before trusting the split."
            )
        for c in comps:
            if c is nuc:
                continue
            if c.peak_depth > nuc.peak_depth and c.span_bp > nuc.span_bp * 0.5:
                notes.append(
                    f"{c.label} is deeper than the nuclear compartment but spans a "
                    f"similar amount of sequence ({c.span_bp / 1000:.0f} kb). That can "
                    f"happen with a second organism in the culture rather than a "
                    f"high copy element."
                )

    if not comps:
        notes.append("No depth humps were found at all. The read set may be too small.")

    return notes
