"""Assembling the small circles, from scratch.

Once triage has pulled the high copy reads out, what is left is a very small
assembly problem with an enormous amount of data: a circle of a few tens of
kilobases covered hundreds or thousands of times over. That is the friendliest
assembly problem there is, and it is small enough that writing the assembler
ourselves is sensible rather than silly.

It is also a problem that general assemblers handle awkwardly. They are built to
produce linear contigs, so a circle usually comes out as a line with its two ends
saying the same thing, and working out whether it really closed is left to you.
Here the molecules genuinely are circles, so closing them is part of the job
rather than an afterthought.

## How the walk knows when to stop

This is the part worth reading, because it is where the first two attempts at
this module went wrong.

The assembler picks a seed read and walks outwards, repeatedly stepping to
whichever overlapping read reaches furthest. On a chromosome that walk ends when
it runs out of reads. On a circle it never ends. Coverage is enormous, so after
going right round the walk arrives back where it started and finds a completely
different set of reads covering the same sequence. It then goes round again. The
first working version of this code happily lapped a twelve kilobase circle more
than thirty times and produced a four hundred kilobase answer.

So the walk watches for the seed. At every step it asks whether the read it is
standing on overlaps the seed read, and whether it has travelled far enough that
this cannot simply be the overlap it set out from. When both are true it has gone
right round, and the distance it has covered is the length of the circle.

That measured length is worth more than inferring circularity afterwards by
comparing the two ends of a draft, because it comes from the overlaps themselves
and it is available before the consensus is built.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass, field

from . import nagcore
from .seqio import reverse_complement


@dataclass
class Placement:
    """One read, laid down at an offset along the draft, on a given strand."""

    read_index: int
    offset: int
    forward: bool
    anchors: int


@dataclass
class CircularAssembly:
    """What came out of assembling one high copy compartment."""

    sequence: str
    is_circular: bool
    n_reads_used: int
    layout: list[Placement] = field(default_factory=list)
    mean_coverage: float = 0.0
    circumference: int = 0
    notes: list[str] = field(default_factory=list)

    def __len__(self) -> int:
        return len(self.sequence)


def downsample(reads: list[str], target_depth: float, expected_bp: int) -> list[str]:
    """Keep roughly enough reads for the depth asked for.

    Longer reads are preferred, because a long read spans more of a repeat and
    the repeats inside the rDNA circle are the hard part of this assembly.

    But not without a limit, and the limit matters. A read longer than the
    element itself has gone right round the circle and started again, so it
    describes the same place twice. The layout here works in straight line
    coordinates, so such a read pulls the walk in two directions at once.

    On the real mitochondrial reads this was not a hypothetical. Asking for
    sixty fold coverage of a fifty kilobase element picked the fifty seven
    longest reads available, which ran from 34 kb to 84 kb with a median of
    38 kb. Three of them were longer than the whole element. The walk took a
    wrong turn and produced a 61 kb circle instead of a 49 kb one, and nothing
    about the result looked wrong from the inside.

    So reads that cannot wrap are used first, and the ones that can are only
    brought in if there are not enough of the others.
    """
    if expected_bp <= 0:
        return list(reads)

    budget = target_depth * expected_bp
    # A little under the element size, so a read cannot quite go round.
    wrap_limit = int(expected_bp * 0.9)
    safe = sorted(
        (r for r in reads if len(r) <= wrap_limit), key=len, reverse=True
    )
    rest = sorted((r for r in reads if len(r) > wrap_limit), key=len, reverse=True)

    kept: list[str] = []
    total = 0
    for read in safe + rest:
        if total >= budget:
            break
        kept.append(read)
        total += len(read)
    return kept


def _place_neighbour(
    anchor_offset: int,
    anchor_forward: bool,
    anchor_len: int,
    other_len: int,
    same_strand: bool,
    shift: int,
    k: int,
) -> tuple[int, bool]:
    """Where a read sits once its partner has already been placed.

    All the fiddly bookkeeping of this module lives in this one function, so it
    is worth being explicit about what the numbers mean.

    The Rust core reports, for a pair, whether they run the same way and a shift.
    For a same strand pair the shift says that position p in the second read is
    position p + shift in the first. For a flipped pair the two coordinate
    systems run against each other, so position p in the second read is position
    shift - p in the first.

    On top of that the already placed read may itself have been laid down
    backwards, which flips the meaning of everything hanging off it. That gives
    four cases, and all four are below.

    The two flipped cases carry a `k` in them and the two same strand cases do
    not. That is not a fudge factor. A shift relates the *start* positions of two
    k-mers, and when a pair is flipped the k-mer that starts at a given place in
    one read ends at the matching place in the other. The length of the k-mer is
    exactly the difference between those two things. Leaving it out puts every
    reverse complemented read fifteen or so bases out of position, which is
    enough to turn a consensus into noise while leaving the overall length
    looking perfectly reasonable.

    Returns the offset of the new read and whether it is laid down forwards.
    """
    if anchor_forward:
        if same_strand:
            # The second read starts `shift` bases along from the first.
            return anchor_offset + shift, True
        # Flipped, so the second read goes down reverse complemented.
        return anchor_offset + shift - other_len + k, False

    # The anchor is itself reversed, so its original coordinates run backwards
    # along the draft and both cases above turn around.
    if same_strand:
        return anchor_offset + anchor_len - other_len - shift, False
    return anchor_offset + anchor_len - shift - k, True


def _overlap_graph(
    reads: list[str],
    k: int,
    w: int,
    min_entropy: float,
    min_anchors: int,
    min_shared: int,
) -> dict[int, list[tuple[int, bool, int, int]]]:
    """Neighbour lists: read -> [(other, same strand, shift, anchors), ...].

    Each overlap is stored from both ends. Seen from the other side a same
    strand shift runs the other way, so it changes sign. A flipped shift is a
    sum rather than a difference, which makes it its own inverse, so it carries
    over unchanged.
    """
    index = nagcore.OverlapIndex(k=k, w=w, min_entropy=min_entropy)
    index.add_reads(reads)
    hits = index.place_overlaps(
        min_shared=min_shared, max_occurrence=10_000, bin_width=64, min_anchors=min_anchors
    )

    graph: dict[int, list[tuple[int, bool, int, int]]] = {}
    for a, b, same_strand, shift, anchors, _, _, _, _ in hits:
        graph.setdefault(a, []).append((b, same_strand, shift, anchors))
        back = -shift if same_strand else shift
        graph.setdefault(b, []).append((a, same_strand, back, anchors))
    for lst in graph.values():
        lst.sort(key=lambda t: t[3], reverse=True)
    return graph


def _walk(
    seed: int,
    reads: list[str],
    graph: dict[int, list[tuple[int, bool, int, int]]],
    rightwards: bool,
    visited: set[int],
    max_span: int,
    k: int,
) -> tuple[list[Placement], int]:
    """Extend a chain of reads out from a seed, one read at a time.

    At each step we look at every read overlapping the current one and keep the
    one that reaches furthest while still overlapping. Taking the furthest reach,
    rather than the strongest overlap, is what keeps the chain moving instead of
    shuffling sideways through a pile of reads that all start in nearly the same
    place.

    ## Spotting the lap

    On a circle the walk never runs out of reads. Coverage is deep, so after
    going right round it meets a fresh set of reads covering the sequence it has
    already crossed, and it carries on into a second lap.

    The way to notice is to watch for a read we have already put down. Every read
    the walk places keeps its offset. When the current read overlaps one of those
    and the overlap says it should sit a long way back from where we actually put
    it, the walk has come round and the gap between the two answers is the length
    of the circle.

    An earlier version only watched for the seed read itself. That works when the
    seed happens to be covered at the point the walk arrives back, and fails
    quietly when it is not: the walk sails past its own starting point and closes
    on the second lap instead, giving a circle exactly twice the real length. The
    rDNA circle came out at 29.9 kb that way, against a true length near 15.8 kb.
    Watching every placed read rather than just the seed fixes it, because by
    then there are dozens of reads laid down near the origin and the walk only
    has to meet one of them.

    Returns the chain and the circle length, or zero if it did not close.
    """
    out: list[Placement] = []
    current = seed
    offset = 0
    forward = True
    placed: dict[int, int] = {seed: 0}

    # A lap has to be longer than this to be believed. Anything shorter is
    # ordinary disagreement between two overlaps rather than a trip round.
    #
    # The bar is deliberately low and fixed, and it took two wrong answers to
    # work out that it should be.
    #
    # The temptation is to scale it to the reads, on the grounds that a lap
    # should be longer than a read. Both ways of doing that fail. Scaling to the
    # longest read put the bar above the size of the rDNA circle, because that
    # read set held a few reads of over 80 kb, so the real lap never counted and
    # the walk went round twice. Scaling to the typical read then did the same
    # thing to the mitochondrion, whose reads are long next to it.
    #
    # The bar does not need to be anywhere near a read length. A read that has
    # been placed correctly sits where the overlap says it should, so the gap is
    # a few tens of bases. Only a genuine trip round the circle produces a gap of
    # thousands. A couple of kilobases clears the jitter with room to spare, and
    # leaves every element larger than that findable.
    min_lap = 2_000

    while True:
        for other, same_strand, shift, _anchors in graph.get(current, ()):
            previous = placed.get(other)
            if previous is None:
                continue
            implied, _ = _place_neighbour(
                offset, forward, len(reads[current]), len(reads[other]), same_strand, shift, k
            )
            gap = implied - previous
            if not rightwards:
                gap = -gap
            if gap > min_lap:
                return out, gap

        best: tuple[int, int, bool, int] | None = None
        for other, same_strand, shift, anchors in graph.get(current, ()):
            if other in visited:
                continue
            new_offset, new_forward = _place_neighbour(
                offset, forward, len(reads[current]), len(reads[other]), same_strand, shift, k
            )
            advance = new_offset - offset
            if not rightwards:
                advance = -advance
            # It has to move us along, and it has to still overlap what we have.
            if advance <= 0 or advance >= len(reads[current]):
                continue
            # "Furthest" means furthest in whichever direction we are walking,
            # so going left we compare negated offsets.
            reach = new_offset if rightwards else -new_offset
            best_reach = None if best is None else (best[1] if rightwards else -best[1])
            if best_reach is None or reach > best_reach:
                best = (other, new_offset, new_forward, anchors)

        if best is None:
            return out, 0
        other, new_offset, new_forward, anchors = best
        if abs(new_offset) > max_span:
            return out, 0
        visited.add(other)
        placed[other] = new_offset
        out.append(Placement(other, new_offset, new_forward, anchors))
        current, offset, forward = other, new_offset, new_forward


def _build_layout(
    reads: list[str],
    k: int,
    w: int,
    min_entropy: float,
    min_anchors: int,
    min_shared: int,
    max_span: int,
    n_seeds: int = 11,
    tolerance: float = 0.05,
) -> tuple[list[Placement], int]:
    """Lay the reads out on one coordinate line.

    This happens in two parts, for two different reasons.

    First a backbone is walked out from a seed read. The backbone decides the
    shape and the length of the answer, so it wants the smallest number of reads
    that reach furthest, each joined to the next by a solid overlap. If the walk
    comes back round to its seed then the element is a circle and we know its
    length. If it does not, the walk is repeated leftwards to pick up the other
    side of what must be a linear element.

    Then every remaining read is dropped onto the backbone wherever it overlaps.
    Those reads add nothing to the shape, but they are what the consensus votes
    with, and a vote taken over three hundred reads is worth a great deal more
    than one taken over the dozen that happened to form the chain.

    ## Why the walk is run more than once

    A greedy walk takes the longest step available at each point and never looks
    back, so one bad overlap sends it on a detour it cannot undo. On the real
    mitochondrial reads this happened for some sets of reads and not others. Ten
    settings out of twelve gave a circle of about 49 kb, matching the published
    genome, and two gave 60 kb, because the walk went down a wrong turn and came
    back round the long way.

    Nothing about the 60 kb answer looks wrong from the inside. It closes into a
    circle, the reads agree with each other along it, and the coverage is even.

    Five starting points turned out not to be enough on that data and eleven
    was, so the walk is run from several different starting reads and the
    lengths are compared. Agreement between walks that set out from different places is
    worth much more than any single walk's confidence in itself, and a wrong
    turn is usually particular to one starting point. The answer kept is the one
    the most walks agree on, and ties go to whichever used the most reads.

    Returns the layout and the circle length, which is zero for a line.
    """
    graph = _overlap_graph(reads, k, w, min_entropy, min_anchors, min_shared)
    if not graph:
        return [], 0

    support = {r: sum(n[3] for n in graph[r]) for r in graph}
    seeds = sorted(support, key=lambda r: support[r], reverse=True)[: max(1, n_seeds)]

    attempts: list[tuple[int, list[Placement], set[int]]] = []
    for seed in seeds:
        visited = {seed}
        backbone = [Placement(seed, 0, True, 0)]
        right, circumference = _walk(seed, reads, graph, True, visited, max_span, k)
        backbone += right
        if not circumference:
            left, _ = _walk(seed, reads, graph, False, visited, max_span, k)
            backbone += left
        attempts.append((circumference, backbone, visited))

    best = _agree_on_an_attempt(attempts, tolerance)
    circumference, backbone, _ = best

    placed = {p.read_index: p for p in backbone}

    # Hang every other read off whichever placed read it overlaps best.
    #
    # The shift stored in graph[x] is written from x's point of view, so for an
    # unplaced read it runs the wrong way round: we need the shift as the read
    # we are hanging it off would see it. Inverting is the same rule used when
    # the graph was built, a sign change for a same strand pair and nothing at
    # all for a flipped one, since a flipped shift is a sum and is its own
    # inverse.
    for read_id, neighbours in graph.items():
        if read_id in placed:
            continue
        for other, same_strand, shift, anchors in neighbours:
            anchor = placed.get(other)
            if anchor is None:
                continue
            from_anchor = -shift if same_strand else shift
            offset, forward = _place_neighbour(
                anchor.offset,
                anchor.forward,
                len(reads[other]),
                len(reads[read_id]),
                same_strand,
                from_anchor,
                k,
            )
            placed[read_id] = Placement(read_id, offset, forward, anchors)
            break

    layout = sorted(placed.values(), key=lambda p: p.offset)
    base = layout[0].offset
    for p in layout:
        p.offset -= base
    return layout, circumference


def _agree_on_an_attempt(
    attempts: list[tuple[int, list[Placement], set[int]]], tolerance: float
) -> tuple[int, list[Placement], set[int]]:
    """Pick the walk that the most other walks agree with.

    Two walks agree when their circle lengths are within `tolerance` of each
    other. The biggest group of agreeing walks wins, and within that group we
    keep the walk that used the most reads.

    Walks that did not close are only considered when no walk closed at all. A
    circle that closed is a much stronger statement than a line that simply ran
    out of reads, so one closed walk beats any number of open ones.
    """
    closed = [a for a in attempts if a[0] > 0]
    pool = closed or attempts
    if not pool:
        return 0, [], set()
    if len(pool) == 1:
        return pool[0]

    best: tuple[int, int, tuple] | None = None  # votes, reads, attempt
    for attempt in pool:
        length = attempt[0]
        if length > 0:
            votes = sum(
                1 for other in pool if abs(other[0] - length) <= tolerance * length
            )
        else:
            votes = sum(1 for other in pool if other[0] == 0)
        key = (votes, len(attempt[1]))
        if best is None or key > (best[0], best[1]):
            best = (votes, len(attempt[1]), attempt)
    return best[2]


def _stitch(reads: list[str], layout: list[Placement], circumference: int = 0) -> str:
    """Build the draft by laying backbone reads end to end.

    Each read contributes the stretch that the reads before it did not already
    cover, so the draft is raw read sequence joined at the seams rather than
    anything averaged.

    That means the draft carries the error rate of the reads themselves, and on a
    noisy long read platform that is a few percent. It gets cleaned up afterwards
    by polishing with a real aligner, which is the same division of labour
    miniasm and racon use.

    The alternative, voting base by base across every read at its offset, is
    implemented below as `_consensus` and is better when the reads are accurate.
    It does not work on indel rich reads, because the offsets come from anchor
    chains and two reads of the same stretch genuinely differ in length once
    indels are in play. Votes then land a few bases apart and the majority is
    meaningless. Stitching sidesteps that by never asking two reads to agree on a
    coordinate.
    """
    if not layout:
        return ""
    pieces: list[str] = []
    covered = 0
    for p in sorted(layout, key=lambda q: q.offset):
        seq = reads[p.read_index]
        if not p.forward:
            seq = reverse_complement(seq)
        end = p.offset + len(seq)
        if end <= covered:
            continue
        start_in_read = max(0, covered - p.offset)
        pieces.append(seq[start_in_read:])
        covered = end
        if circumference and covered >= circumference:
            break
    draft = "".join(pieces)
    if circumference:
        draft = draft[:circumference]
    return draft


def _consensus(
    reads: list[str], layout: list[Placement], circumference: int = 0
) -> tuple[str, float]:
    """Vote base by base across the laid out reads.

    Every read contributes one vote at each position it covers. The winning base
    at each position is kept, and positions nobody covers become N.

    When `circumference` is given, the draft is that long and a read running off
    the end wraps back to the start, because on a circle it really does. That
    also means every read votes, including the ones the walk met on its way back
    round to the seed.

    This is a simple vote and it is deliberately simple. The offsets come from
    anchor chains, so they can be out by an indel here and there, and a cleverer
    scheme built on slightly wrong offsets would not be more correct, only more
    confident. Polishing with a real aligner is the right way to sharpen this,
    and the pipeline does that afterwards.
    """
    if not layout:
        return "", 0.0

    length = circumference or max(p.offset + len(reads[p.read_index]) for p in layout)
    if length <= 0:
        return "", 0.0

    votes: list[Counter] = [Counter() for _ in range(length)]
    for p in layout:
        seq = reads[p.read_index]
        if not p.forward:
            seq = reverse_complement(seq)
        start = p.offset
        for i, base in enumerate(seq):
            if base not in "ACGTacgt":
                continue
            pos = start + i
            if circumference:
                pos %= circumference
            elif pos >= length:
                break
            votes[pos][base.upper()] += 1

    out: list[str] = []
    covered = 0
    depth_total = 0
    for column in votes:
        if column:
            out.append(column.most_common(1)[0][0])
            covered += 1
            depth_total += sum(column.values())
        else:
            out.append("N")
    mean_cov = depth_total / covered if covered else 0.0
    return "".join(out), mean_cov


def _trim_end_repeat(seq: str, k: int, min_overlap: int) -> int:
    """How many bases the end of a draft repeats from its start.

    Only used for drafts the walk could not close. A linear draft of a circular
    molecule says the same thing at both ends, so finding that agreement is a
    second, weaker way to spot a circle. The search runs from the longest
    plausible repeat downwards and takes the first clean match, since a long
    agreement is much harder to get by chance than a short one.
    """
    n = len(seq)
    top = n // 2
    for length in range(top, max(min_overlap, k) - 1, -1):
        if seq[:length] == seq[n - length :]:
            return length
    return 0


def rotate_to_canonical_start(seq: str, k: int = 16) -> str:
    """Spin a circle round so it always starts at the same place.

    A circle has no beginning, so the same molecule assembled twice can come out
    starting anywhere, and two files that look nothing alike can describe exactly
    the same thing. We pick a starting point from the sequence itself: the
    position of the smallest k-mer, which does not depend on where the assembly
    happened to begin. Ties go to the earliest position, so the answer is always
    the same.
    """
    if len(seq) <= k:
        return seq
    doubled = seq + seq[: k - 1]
    best_pos = 0
    best = doubled[:k]
    for i in range(1, len(seq)):
        candidate = doubled[i : i + k]
        if candidate < best:
            best = candidate
            best_pos = i
    return seq[best_pos:] + seq[:best_pos]


def assemble_circle(
    reads: list[str],
    expected_bp: int = 0,
    target_depth: float = 60.0,
    k: int = 15,
    w: int = 10,
    min_entropy: float = 1.0,
    min_anchors: int = 4,
    min_shared: int = 3,
    rotate: bool = True,
    vote: bool | None = None,
    n_seeds: int = 11,
) -> CircularAssembly:
    """Assemble one high copy compartment and try to close it into a circle.

    `expected_bp` is the span triage estimated for this compartment. It decides
    how many reads to keep and how far the walk is allowed to run before it is
    treated as having gone wrong. Zero means use every read and set no limit.

    `vote` picks how the draft sequence is built. Voting base by base across all
    the reads gives a clean answer when the reads are accurate, which is the case
    for short reads and for HiFi. On noisy long reads it does not work, because
    indels mean two reads of the same stretch disagree about how long it is, so
    the votes land out of register. The default follows that: vote when the reads
    look accurate, otherwise stitch the backbone together and leave the cleaning
    up to the polishing step.
    """
    notes: list[str] = []
    reads = [r for r in reads if len(r) >= k * 2]
    if not reads:
        return CircularAssembly("", False, 0, notes=["no reads long enough to assemble"])

    kept = downsample(reads, target_depth, expected_bp) if expected_bp else list(reads)
    if len(kept) < len(reads):
        notes.append(
            f"used {len(kept)} of {len(reads)} reads, which is about "
            f"{target_depth:.0f} fold coverage of an element of {expected_bp:,} bp"
        )

    # Overlaps are measured on the reads themselves, not on a collapsed copy of
    # them. Collapsing homopolymers finds pairs well, but the offsets it reports
    # are in collapsed coordinates, and scaling those back is only ever
    # approximate. That is no use when the offsets are what the consensus votes
    # on.
    max_span = max(expected_bp * 4, 200_000) if expected_bp else 2_000_000
    layout, circumference = _build_layout(
        kept, k, w, min_entropy, min_anchors, min_shared, max_span, n_seeds=n_seeds
    )

    if not layout:
        return CircularAssembly(
            "", False, 0, notes=notes + ["no overlaps found between these reads"]
        )
    if len(layout) < len(kept):
        notes.append(
            f"{len(kept) - len(layout)} of {len(kept)} reads could not be joined to "
            f"the rest and were left out"
        )

    if vote is None:
        # Short reads come off accurate platforms. Long ones usually do not.
        vote = sum(len(r) for r in kept) / len(kept) < 1_000

    if vote:
        draft, mean_cov = _consensus(kept, layout, circumference)
    else:
        draft = _stitch(kept, layout, circumference)
        total = sum(len(kept[p.read_index]) for p in layout)
        mean_cov = total / len(draft) if draft else 0.0
        notes.append(
            "the draft is raw read sequence stitched together, so it still carries "
            "the error rate of the reads. Polish it before using it"
        )
    if not draft:
        return CircularAssembly("", False, 0, notes=notes + ["consensus came out empty"])

    is_circular = circumference > 0
    if is_circular:
        sequence = draft
        notes.append(
            f"the walk came back round to its starting read after {circumference:,} bp, "
            f"so this element is a circle"
        )
    else:
        # The walk did not close. Fall back to looking for the two ends of the
        # draft agreeing, which is the older and weaker test.
        repeat = _trim_end_repeat(draft, k, min_overlap=max(k, 40))
        if repeat:
            sequence = draft[: len(draft) - repeat]
            circumference = len(sequence)
            is_circular = True
            notes.append(
                f"the walk did not close, but the two ends of the draft repeat "
                f"{repeat:,} bp of each other, so this is reported as a circle"
            )
        else:
            sequence = draft
            notes.append(
                "the two ends do not meet, so this is reported as a line. That can "
                "mean the element is not circular, or that coverage does not reach "
                "all the way round"
            )

    if is_circular and rotate:
        sequence = rotate_to_canonical_start(sequence)

    n_gaps = sequence.count("N")
    if n_gaps:
        notes.append(f"{n_gaps} positions had no read covering them and are N")

    return CircularAssembly(
        sequence=sequence,
        is_circular=is_circular,
        n_reads_used=len(layout),
        layout=layout,
        mean_coverage=mean_cov,
        circumference=circumference,
        notes=notes,
    )


def repeat_content(seq: str, k: int = 50, min_copies: int = 2) -> float:
    """Share of a sequence sitting inside a k-mer that occurs more than once.

    The rDNA circle is a good deal more repetitive than most small genomes, and
    that is worth reporting next to the assembly, because repeats are what limits
    how much of it short reads can resolve.
    """
    if len(seq) < k:
        return 0.0
    counts: Counter = Counter()
    windows = len(seq) - k + 1
    for i in range(windows):
        counts[seq[i : i + k]] += 1
    repeated = sum(1 for i in range(windows) if counts[seq[i : i + k]] >= min_copies)
    return repeated / windows


def resolvability_note(seq: str, reads: list[str], k: int = 50) -> str:
    """A plain sentence on whether the reads can get across this element's repeats.

    A repeat can only be resolved by a read that spans it with unique sequence on
    both sides. Short reads cannot do that for the long direct repeats published
    in the Naegleria rDNA circle, so saying so is more useful than quietly
    producing a confident answer.
    """
    share = repeat_content(seq, k=k)
    longest = max((len(r) for r in reads), default=0)
    if share < 0.05:
        return (
            f"About {share:.0%} of this element is repeated, so repeats are not a "
            f"problem here."
        )
    if longest >= len(seq):
        return (
            f"About {share:.0%} of this element is repeated, but the longest read "
            f"({longest:,} bp) spans the whole element, so the repeats are resolved."
        )
    return (
        f"About {share:.0%} of this element sits inside a repeat, and the longest read "
        f"is {longest:,} bp against an element of {len(seq):,} bp. Any repeat longer "
        f"than a read cannot be resolved from this data, so treat the copy number and "
        f"the order of the repeated blocks as unsettled."
    )


def reads_from_element(
    reads: list[str],
    sequence: str,
    is_circular: bool = True,
    k: int = 15,
    w: int = 10,
    min_share: float = 0.05,
) -> list[bool]:
    """Which reads came from an element we have already assembled.

    Each read is sketched down to its anchor k-mers and we ask what share of them
    the element also has. A read from the element shares a good fraction even
    when both are noisy. A read from somewhere else shares essentially none, so
    the two groups separate cleanly and the threshold does not need to be
    delicate.

    On real nanopore reads against the mitochondrial genome this came out as a
    quarter of reads at a share of zero and the rest spread between 0.05 and
    0.46, which is about as clean a split as this kind of test gives.

    For a circular element the sequence is doubled first, so that anchors lying
    across the join are not missed.
    """
    target = sequence + sequence if is_circular else sequence
    present = {h for h, _, _ in nagcore.sketch(target, k, w, 1.0)}
    out: list[bool] = []
    for read in reads:
        anchors = nagcore.sketch(read, k, w, 1.0)
        if not anchors:
            out.append(False)
            continue
        share = sum(1 for h, _, _ in anchors if h in present) / len(anchors)
        out.append(share >= min_share)
    return out


def same_element(a: str, b: str, k: int = 15, w: int = 10, min_share: float = 0.4) -> bool:
    """Whether two assembled sequences describe the same molecule.

    Compared through their anchor k-mers, so two drafts of one element match
    even though they start at different points, run in different directions and
    carry different sequencing errors.

    The share is measured against the shorter of the two, so a fragment of an
    element counts as the same thing rather than as something new.
    """
    if not a or not b:
        return False
    anchors_a = {h for h, _, _ in nagcore.sketch(a + a, k, w, 1.0)}
    anchors_b = {h for h, _, _ in nagcore.sketch(b + b, k, w, 1.0)}
    if not anchors_a or not anchors_b:
        return False
    shared = len(anchors_a & anchors_b)
    return shared / min(len(anchors_a), len(anchors_b)) >= min_share


def assemble_elements(
    reads: list[str],
    expected_bp: int = 0,
    max_elements: int = 4,
    min_reads: int = 20,
    min_length_bp: int = 1_000,
    **kwargs,
) -> list[CircularAssembly]:
    """Assemble every element in a compartment, not just the first one.

    One compartment can hold more than one molecule. Naegleria is the obvious
    case: the mitochondrial genome and the rDNA circle are both present in many
    copies, and depending on the library their depths can land close enough
    together that the triage step sees one hump rather than two.

    So after each element is assembled, the reads that came from it are taken out
    and whatever is left goes round again. This stops when nothing more can be
    assembled, when the elements coming out are too small to be real, or when
    `max_elements` is reached.

    Returns the elements in the order they were assembled, which is deepest and
    best covered first.
    """
    remaining = list(reads)
    found: list[CircularAssembly] = []

    for _ in range(max_elements):
        if len(remaining) < min_reads:
            break
        assembly = assemble_circle(remaining, expected_bp=expected_bp, **kwargs)
        if not assembly.sequence or len(assembly.sequence) < min_length_bp:
            break

        # The reads belonging to an element are not always all claimed, so the
        # leftovers can assemble into the same thing a second time. Say so and
        # stop rather than reporting one molecule twice.
        if any(same_element(assembly.sequence, seen.sequence) for seen in found):
            break

        found.append(assembly)
        belongs = reads_from_element(
            remaining, assembly.sequence, is_circular=assembly.is_circular
        )
        left = [r for r, used in zip(remaining, belongs) if not used]
        assembly.notes.append(
            f"{len(remaining) - len(left):,} of {len(remaining):,} reads came from "
            f"this element and were set aside"
        )
        if len(left) == len(remaining):
            # Nothing was claimed, so going round again would only find the same
            # thing a second time.
            break
        remaining = left
        # Later elements are different sizes, so the earlier estimate no longer
        # applies and we let the assembler use everything it has.
        expected_bp = 0

    return found
