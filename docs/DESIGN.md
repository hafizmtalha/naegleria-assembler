# Why this assembler exists

Most assemblers are built with a common kind of genome in mind. A bacterium is one
circle. A human is a big diploid genome with a lot of repeats. A virus is tiny. The
tools work well on those shapes because the shape is baked into their assumptions.

Naegleria is a different shape, and the difference is not small.

## The thing that breaks normal assemblers

Naegleria keeps its ribosomal RNA genes off the chromosomes. They sit on their own
small circle of DNA, about 15.8 kb, and the cell holds roughly four thousand copies of
it. The mitochondrial genome is a second circle, about 49 kb, and it is also present in
many copies.

Now think about what that does to a pile of sequencing reads. The nuclear genome is
about 29.5 Mb and sits at some baseline depth. The small circles sit at hundreds or
thousands of times that depth. They are a rounding error in total length and a huge
share of the reads.

This is not a guess. The team that produced the 2019 nanopore assembly ran into exactly
this. They wrote that these high coverage elements made the combined assembly perform
poorly, so they pulled those reads out and assembled them on their own.

Their way of pulling them out was to map reads against a known rDNA circle from a
related species and throw away whatever stuck. That works, but it needs you to already
have the answer. If you are assembling a new isolate, or a Naegleria relative nobody has
sequenced, you do not have a reference circle to map against.

So the first thing this tool does is find those compartments without a reference.

## Idea one: sort the reads by depth before assembling anything

Every read gets a depth score. We count canonical k-mers across the whole read set,
then give each read the median count of its own k-mers. Median, not mean, so one
strange region does not drag a read around.

Plot those scores and you do not get one hump. You get several. The nuclear reads make
a broad hump at the baseline. The mitochondrial reads make a hump far to the right. The
rDNA circle makes a hump further right still. Anything from a bacterium in the culture
sits off on its own, usually at a different depth and a different GC.

We fit those humps and hand each read a label. Nuclear, high copy circle, or odd one
out. Nothing in that needs a reference genome. It only needs the reads.

GC comes in as a second axis, because depth alone can be ambiguous. A contaminant at
the same depth as the mitochondrion is easy to separate once you notice its base
composition is nothing like Naegleria.

## Idea two: assemble each compartment the way that compartment deserves

Once the reads are sorted, the problem splits into pieces that are each much easier
than the whole.

The small circles are now tiny problems with enormous coverage. That is the friendliest
assembly problem there is, and it is small enough that writing our own assembler for it
is sensible rather than silly. We build an overlap graph, find the cycle, and rotate the
answer to a fixed starting point so two runs give the same string.

Circularity matters here and normal assemblers often leave it to you to work out. These
elements really are circles, so the tool checks for the overlap between the two ends and
reports whether it closed.

The rDNA circle has a complication worth naming. A little under half of its non
ribosomal part is made of repeated blocks, some over two kilobases long. Short reads
cannot get across those. The tool says so plainly rather than quietly producing a
confident wrong answer, and in short read only mode it reports that part as unresolved.

The nuclear reads, with the high copy noise taken out, go to the nuclear assembler. In
stage one that is a proven tool. In stage two it is ours.

## Idea three: stop letting poly-A ruin the anchors

Naegleria is AT rich. GC sits around a third of the genome.

Assemblers pick anchor k-mers to find overlaps between reads, and the usual way to pick
them takes the smallest k-mer in a window by alphabetical order. A runs before C, G and
T. In an AT rich genome that quietly means you keep choosing runs of A, which are
everywhere, carry almost no information about where you are, and are the single worst
thing to try to line two reads up on.

We pick anchors by a hash instead of by the alphabet, and we throw out any anchor whose
sequence is too low in complexity to be worth trusting. On top of that, long read runs
through homopolymers are the classic error, so overlaps are found on a homopolymer
compressed view of the read while the original sequence is kept for the final consensus.

## What runs in what language

The hot loops are in Rust. Counting k-mers over a few gigabytes of reads, sketching
them, and comparing sketches are the parts where Python would be too slow to be useful.

Everything above that is Python, because the point of this repo is that someone can read
it, follow the logic and change it for their own organism.

## Staging

Stage one is the compartment work plus our own circular assembler, with a proven tool
doing the nuclear assembly. Stage two swaps in our own nuclear assembler as a second
engine you can select.

The honest reason for that order is that the compartment triage is where the real gain
is. A from scratch nuclear assembler is a good thing to have read and understood, but it
is not going to beat a mature tool, and pretending otherwise would be silly.
