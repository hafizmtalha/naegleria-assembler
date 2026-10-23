# nagasm

A genome assembler built around the way Naegleria keeps its DNA.

Most assemblers assume that coverage is roughly even across a genome. Naegleria
breaks that assumption, and it breaks it badly enough that general tools
struggle. This one sorts reads into compartments by depth before it assembles
anything, so the small circles that Naegleria carries in thousands of copies stop
getting in the way of the chromosomes.

Everything here runs on public data and you can reproduce the results below with
two commands.

![Depth spectrum](docs/figures/real_depth_spectrum.png)

That plot is the whole idea in one picture. The big hump on the left is the
nuclear genome. The spike on the right is the mitochondrial genome, sitting at
seventy times the depth. On a normal linear plot you would never see it.

## The problem

Naegleria does not keep its ribosomal RNA genes on its chromosomes. They sit on
their own small circle of DNA, and the cell holds something like four thousand
copies of it. The mitochondrial genome is a second circle, also in many copies.

Think about what that does to a pile of reads. The nuclear genome is about 29.5
Mb and sits at whatever depth you sequenced to. The circles are a fraction of a
percent of that length, but because of their copy number they can account for a
large share of the reads. In the public dataset used below, one circle of 49 kb
pulled in more than one read in six.

Hand that to an ordinary assembler and it has a hard time, because almost
everything an assembler does assumes that unusually deep sequence is a repeat to
be collapsed.

This is not a guess. The team who produced the 2019 nanopore assembly of this
organism ran into exactly this. They wrote that these high coverage elements made
the combined assembly perform poorly, so they pulled those reads out and
assembled them separately.

Their way of pulling them out was to map reads against a known rDNA circle from a
related species and throw away whatever stuck. That works, but it needs you to
already have the answer. If you are working on a new isolate, or a relative
nobody has sequenced yet, there is no reference circle to map against.

So this tool finds the compartments in the reads themselves.

## What is new here

### 1. Sorting reads by depth, with no reference

Count k-mers across the whole read set, then give every read the median count of
its own k-mers. Reads from a four thousand copy circle are built from k-mers seen
far more often than k-mers from a single copy chromosome, so that one number
already separates them.

To decide where the boundaries go we look at the k-mer count histogram. Each
compartment makes its own hump. The useful part is that the area under a hump
tells you how much distinct sequence sits at that depth, which is what tells a
large genome apart from a small element held in many copies.

That area is also what decides which compartment is the nuclear one. It is
simply the compartment covering the most distinct sequence. It is not the tallest
hump, and it is not the one with the most reads. On this organism those are
different things, and picking either of the easy ones gets it wrong.

One detail that turned out to matter a lot. The peak finding runs on the log of
the histogram. The compartments differ in size by three orders of magnitude, so
any rule of the form "ignore peaks shorter than some share of the tallest" throws
away exactly the compartment this tool exists to find. There is a test for that.

### 2. Each compartment assembled the way it deserves

Once the reads are sorted, the problem splits into pieces that are each much
easier than the whole.

The small circles become tiny problems with enormous coverage, which is friendly
enough that writing our own assembler for them is sensible. It finds overlaps,
walks them into a layout, works out whether the walk came back to where it
started, and rotates the answer to a fixed starting point so that two runs agree.

Circularity is measured rather than guessed. The walk watches for the read it set
out from, and the distance it covered getting back there is the length of the
circle. A general assembler usually hands you a line and leaves you to work out
whether its two ends mean it closed.

A compartment can hold more than one molecule, so after each element is
assembled, the reads belonging to it are set aside and whatever is left goes
round again.

The nuclear reads go to Flye or SPAdes. That is on purpose. Those tools are many
years of work on the hard parts of assembly and a hand written replacement would
be worse. There is one of our own in `olc.py` for reading, and you can switch to
it with `--engine nagasm`, but it is there to be understood rather than to
compete.

### 3. Anchors that are not ruined by an AT rich genome

Assemblers pick anchor k-mers to find overlaps between reads. The usual way takes
the smallest k-mer in a window by alphabetical order. A comes before C, G and T.

In a genome that is about two thirds A and T, that quietly means you keep
choosing runs of A. Those are everywhere, they say almost nothing about where you
are, and lining two reads up on one is close to meaningless.

So anchors are picked by a hash instead of by the alphabet, and any k-mer whose
own base composition is too lopsided is thrown out before it is ever considered.
There is a test that poly-A contributes no anchors at all.

## Does it work

### On real public data

The reads are from [PRJNA541227](https://www.ebi.ac.uk/ena/browser/view/PRJNA541227),
the nanopore run behind the 2019 assembly of Naegleria fowleri ATCC 30863. The run
below used the first 94,572 reads, which is 0.63 Gb, or about 21 times over the
nuclear genome. These are older nanopore reads with an error rate around ten
percent, which makes this a hard test rather than an easy one.

**The mitochondrial genome came out whole.**

| | Our assembly | Published |
|---|---|---|
| Length | 48,281 bp | 49,531 bp |
| Shape | closed circle | circle |
| Covers the reference | 97.8% | |
| Identity where aligned | 98.0% | |

That is against the published mitochondrial genome
([NC_021104.1](https://www.ncbi.nlm.nih.gov/nuccore/NC_021104.1)). Nothing
substantial is missing. The two percent that does not match is leftover nanopore
error, and it would come down with more polishing rounds or a handful of short
reads.

Of the reads that really do come from the mitochondrion, 97.1 percent were put in
the right compartment, with no reference used at any point.

**The rDNA circle came out too, and the story behind it is worth telling.**

| | Our assembly | Published |
|---|---|---|
| Length | 15,644 bp | 15,786 bp (strain LEE) |
| Shape | closed circle | circle |
| Covers the reference | 95.7% | |
| Identity where aligned | 98.1% | |

The published circle used for that comparison
([MT741533.1](https://www.ncbi.nlm.nih.gov/nuccore/MT741533.1)) is from a
different strain, so some of the difference is real biology rather than error.

What is interesting is that the depth sorting on its own did not find it. In this
library the rDNA sits at about eleven times the nuclear depth, which sounds like
plenty, but it is only fifteen kilobases of sequence. The hump it makes is tiny,
and it lands in the long tail of nuclear repeats that runs up from the nuclear
peak. It is a shoulder, not a peak.

Lowering the bar for what counts as a peak does not fix this. It shatters the
nuclear peak into a dozen imaginary compartments long before it finds the real
one. We tried, and it is in the git history.

What works is giving up on finding a shape and instead collecting every read more
than a few times deeper than the nuclear peak, then letting the assembler sort it
out. Most of what comes back is nuclear repeat and assembles into unremarkable
straight pieces. The one thing that closed into a circle was the rDNA.

That step only exists because the real data refused to behave like the synthetic
data. It is in `TriageResult.deep_outlier_indices`, with the reasoning written
next to it.

**Three smaller circles also came out, and we do not know what they are.**

They are about 3.3 kb each, they close cleanly, they are covered deeply, and
they share no sequence at all with either of the two circles above or with the
published mitochondrial genome. They may be real. They may be the walk closing a
loop that is not there. The pipeline reports them with their repeat share next to
them and does not claim they are anything. Saying so is better than quietly
dropping them or quietly promoting them.

**The nuclear assembly stands up next to the published one.**

| | Our assembly | Published |
|---|---|---|
| Total length | 28,093,592 bp | 29.5 Mb |
| Contigs | 65 | 81 to 90 |
| N50 | 714,291 bp | 717,491 bp |

The published assembly used far more data than the 0.63 Gb subset here, so this
is not a claim that our nuclear assembly is better. The point is that the reads
left over after the circles are taken out still build into chromosome sized
pieces, which is what the sorting step is meant to protect.

Checking that assembly properly turned up a real bug in our own tool. A 48,897 bp
contig in the first nuclear run was a complete second copy of the mitochondrion.
The depth sorting had put 97.1 percent of the mitochondrial reads in the right
place, but the three percent it missed still came to roughly sixty times coverage
of a 49 kb molecule, which is plenty for Flye to rebuild it. Being nearly right
is not good enough when the leftovers are that deep.

The fix is to stop trusting depth alone. Once an element is assembled, every read
heading for the nuclear set is checked against it directly, and the ones that
match are dropped. That removed 1,980 reads. In the run above there is no longer
any contig matching the published mitochondrial genome, and the earlier output is
kept in the git history so the difference can be seen.

### On data where we know the answer

Real data tells you whether something works. Synthetic data tells you how well,
because you can compare against the truth.

The test below builds a read set shaped like a Naegleria library, with a
chromosome at ordinary depth and two circles far above it, then runs the whole
pipeline.

| Element | True length | Assembled | Identity |
|---|---|---|---|
| Nuclear | 300,000 bp | 299,986 bp | 99.89% |
| Mitochondrion | 12,000 bp | 12,002 bp | 99.75% |
| rDNA circle | 5,000 bp | 5,021 bp | 99.40% |

Those reads carry a two percent error rate. Every read was put in the right
compartment.

Separation holds up from clean reads to about five percent error, which covers
modern nanopore and PacBio HiFi. Above that it degrades, and the tool says so in
its own output rather than quietly handing you a wrong answer.

## Running it

You need [conda](https://docs.conda.io) or [mamba](https://mamba.readthedocs.io).
A Rust compiler is installed into the environment, so you do not need one
already.

```bash
git clone https://github.com/hafizmtalha/naegleria-assembler.git
cd naegleria-assembler

mamba env create -f environment.yml
conda activate nagasm
maturin develop --release
```

Look at your reads before assembling them. This takes about ten seconds and tells
you whether the compartments are really there.

```bash
nagasm triage reads.fastq.gz --scale 20
```

```
compartment       depth  vs nuclear         span     reads     GC
-----------------------------------------------------------------
nuclear             18x          1x  16,900,240    69,763  42.7%
high_copy_1       1328x         74x      37,460    23,754  33.3%
```

Then run the whole thing.

```bash
nagasm assemble --long reads.fastq.gz -o results --threads 8 --scale 20
```

Short reads instead, or as well:

```bash
nagasm assemble --short r1.fastq.gz r2.fastq.gz -o results
nagasm assemble --long ont.fastq.gz --short r1.fastq.gz r2.fastq.gz -o results
```

There is a Snakemake workflow too, if you prefer one.

```bash
snakemake --cores 8
```

## What you get

```
results/
├── report.md                     what happened, in plain sentences
├── summary.json                  the same thing, for scripting
├── assembly/
│   ├── high_copy_elements.fasta  the circles, each rotated to a fixed start
│   └── nuclear.fasta             the chromosomes
└── figures/
    ├── depth_spectrum.png        the plot at the top of this page
    ├── reads_depth_vs_gc.png     every read as a dot
    └── read_lengths.png          read length by compartment
```

## Settings worth knowing

Everything lives in `config/config.yaml`, and these three are the ones that
actually change the answer.

| Setting | Default | What it does |
|---|---|---|
| `scale` | 100 | Keep about one k-mer in this many. This is the memory knob. **Use 20 for noisy long reads.** See below. |
| `k` | 21 | k-mer length for depth sorting. Drop it if your reads are very noisy. |
| `deep_outlier_factor` | 4 | How many times deeper than the nuclear peak a read has to be before it is looked at on its own. This is what finds the rDNA. |

The `scale` default is wrong for noisy reads and it is worth saying why, because
it cost an evening to work out. At the default of 100, each read contributes only
a few dozen sampled k-mers. If the reads are noisy, most of those carry an error
and get thrown away, which leaves too few to take a median of. On the real
nanopore data the default produced nonsense, and the tool said so in its warnings.
Dropping to 20 fixed it completely.

## Where it falls short

Said plainly, because a tool that hides this is worse than useless.

* **Very noisy reads break the depth sorting.** Above about ten percent error the
  compartments stop separating cleanly. The run prints warnings when this
  happens rather than carrying on quietly.
* **Repeats inside the rDNA circle are not resolved.** Nearly half of the
  non ribosomal part of that circle is made of repeated blocks, some over two
  kilobases long. A repeat can only be settled by a read that spans it with
  unique sequence on both sides. The tool reports whether your reads are long
  enough and says the copy number is unsettled when they are not.
* **Our own nuclear assembler is not competitive.** It is readable and it is
  correct on clean data, but Flye will beat it on anything real. That is why Flye
  is the default.
* **The draft from our circular assembler carries the reads' error rate.** It
  stitches read sequence rather than averaging it, and polishing afterwards is
  what makes it accurate. This is the same split miniasm and racon use.
* **Short reads cannot do this on their own.** They can find the compartments,
  but they cannot get across the repeats inside the circles, so the circles will
  come out broken.

## How it is built

The hot loops are Rust, because counting k-mers across gigabytes of reads and
matching their sketches is not something Python should be asked to do. Everything
that decides what any of it means is Python, so you can read it and change it.

```
src/nagasm/
├── triage.py      sorting reads into compartments. The heart of it
├── circular.py    our own assembler for the small circles
├── olc.py         our own assembler for chromosomes, for reading
├── nuclear.py     handing the nuclear reads to Flye or SPAdes
├── polish.py      minimap2 and racon, including the trick for circles
├── stats.py       N50 and friends
├── report.py      figures and the written report
└── pipeline.py    the order the steps happen in

rust/src/
├── kmer.rs        encoding, canonical k-mers, entropy, homopolymers
├── depth.rs       counting k-mers and describing reads
└── sketch.rs      anchors, and where two reads line up
```

There are 79 Python tests and 18 Rust tests.

```bash
pytest
cd rust && cargo test --no-default-features
```

Quite a few of those tests exist because the code failed them while it was being
written. Three of the bugs produced assemblies that looked entirely reasonable in
the summary numbers and were noise underneath, so those have tests that check the
sequence against a known answer rather than just its length.

## Reading the numbers

N50 rewards long contigs and knows nothing about whether they are right. An
assembler that glues two chromosomes together scores better than one that leaves
them apart. On this organism in particular, an assembly that quietly collapsed
the thousands of copies of the rDNA circle into one long array would look
excellent by that measure and be badly wrong.

Read the contig numbers next to the compartment table, not instead of it.

## References

* Liechti N, Schürch N, Bruggmann R, Wittwer M. Nanopore sequencing improves the
  draft genome of the human pathogenic amoeba Naegleria fowleri. *Scientific
  Reports* 9, 16040 (2019).
* Fritz-Laylin LK et al. The genome of Naegleria gruberi illuminates early
  eukaryotic versatility. *Cell* 140, 631 to 642 (2010).
* Ali IKM et al. Draft chromosome sequences of a clinical isolate of the
  free living ameba Naegleria fowleri. *Microbiology Resource Announcements* 10,
  e01034-20 (2021).
* Kolmogorov M et al. Assembly of long, error prone reads using repeat graphs.
  *Nature Biotechnology* 37, 540 to 546 (2019).
* Vaser R et al. Fast and accurate de novo genome assembly from long uncorrected
  reads. *Genome Research* 27, 737 to 746 (2017).
* Li H. Minimap2: pairwise alignment for nucleotide sequences. *Bioinformatics*
  34, 3094 to 3100 (2018).

## License

MIT. See [LICENSE](LICENSE).
