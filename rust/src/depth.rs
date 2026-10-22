//! How deeply each read is covered, worked out from the reads alone.
//!
//! This is the measurement the whole compartment idea rests on. A read that came
//! from the rDNA circle is made of k-mers the cell holds thousands of copies of,
//! so those k-mers turn up again and again across the read set. A read from a
//! single copy stretch of a chromosome is made of k-mers that turn up far less
//! often. Counting k-mers therefore tells us which compartment a read came from
//! without ever looking at a reference genome.
//!
//! ## Why we only keep a sample of the k-mers
//!
//! Counting every k-mer in a few gigabytes of noisy long reads is not possible in
//! a sane amount of memory. Sequencing errors are the problem. Each error creates
//! a whole run of k-mers that exist nowhere else, so the table fills up with
//! things seen exactly once.
//!
//! So we keep a fixed fraction instead, chosen by hash. A k-mer is kept when its
//! hash falls below a cutoff, which works out to roughly one in `scale` of them.
//! The choice is made from the k-mer itself, so the same k-mer is either always
//! kept or never kept, no matter which read it turns up in. That is what makes
//! the counts comparable across reads.
//!
//! At the default scale of 100 a six kilobase read still contributes around sixty
//! sampled k-mers, which is plenty to take a median of, and the table shrinks by
//! about a hundred fold.

use crate::kmer::{each_canonical_kmer, gc_fraction, splitmix64};
use pyo3::prelude::*;
use rayon::prelude::*;
use rustc_hash::FxHashMap;

/// Counts of the sampled k-mers, plus the settings used to build them.
#[pyclass]
pub struct DepthIndex {
    k: usize,
    scale: u64,
    cutoff: u64,
    counts: FxHashMap<u64, u32>,
}

/// Decide whether a k-mer is in the sample, and give back its hash.
#[inline]
fn sampled(code: u64, cutoff: u64) -> Option<u64> {
    let h = splitmix64(code);
    if h < cutoff {
        Some(h)
    } else {
        None
    }
}

/// Median of a slice, taken by sorting a copy. The slices here are one read's
/// worth of k-mers, so this is cheap and the clarity is worth more than the
/// speed of a selection algorithm.
fn median(values: &mut [u32]) -> f64 {
    if values.is_empty() {
        return 0.0;
    }
    values.sort_unstable();
    let n = values.len();
    if n % 2 == 1 {
        values[n / 2] as f64
    } else {
        (values[n / 2 - 1] as f64 + values[n / 2] as f64) / 2.0
    }
}

impl DepthIndex {
    /// Counting, with no Python involved. `add_reads` wraps this.
    pub fn count_into(&mut self, reads: &[String]) {
        let k = self.k;
        let cutoff = self.cutoff;
        let partial = reads
            .par_chunks(2048)
            .map(|chunk| {
                let mut local: FxHashMap<u64, u32> = FxHashMap::default();
                for read in chunk {
                    each_canonical_kmer(read.as_bytes(), k, |code| {
                        if let Some(h) = sampled(code, cutoff) {
                            *local.entry(h).or_insert(0) += 1;
                        }
                    });
                }
                local
            })
            .reduce(FxHashMap::default, |mut a, mut b| {
                if a.len() < b.len() {
                    std::mem::swap(&mut a, &mut b);
                }
                for (key, value) in b {
                    *a.entry(key).or_insert(0) += value;
                }
                a
            });

        if self.counts.is_empty() {
            self.counts = partial;
        } else {
            for (key, value) in partial {
                *self.counts.entry(key).or_insert(0) += value;
            }
        }
    }

    /// Per read description, with no Python involved. `profile_reads` wraps this.
    ///
    /// `min_count` is the thing that makes this work on noisy reads. A sequencing
    /// error creates k-mers that exist nowhere else, so they sit at a count of one
    /// or two. On a long read platform at ten percent error with k of 21, only
    /// about a ninth of a read's k-mers survive error free, so nearly nine in ten
    /// of them are junk sitting at a count of one. A plain median over all of them
    /// therefore reports one, no matter which compartment the read came from, and
    /// the whole separation collapses.
    ///
    /// Dropping everything below `min_count` first leaves the k-mers that really
    /// are shared across the read set, and the median of those says what the
    /// read's depth actually is. Pass 0 or 1 to keep every k-mer.
    ///
    /// Returns, per read: median depth, mean depth, GC share, length, how many
    /// k-mers were sampled, and how many of those survived the cut. That last
    /// number is what tells you whether to believe the first one.
    pub fn profile_into(
        &self,
        reads: &[String],
        min_count: u32,
    ) -> Vec<(f64, f64, f64, usize, usize, usize)> {
        let k = self.k;
        let cutoff = self.cutoff;
        let counts = &self.counts;
        let floor = min_count.max(1);
        reads
            .par_iter()
            .map(|read| {
                let bytes = read.as_bytes();
                let mut sampled_n = 0usize;
                let mut kept: Vec<u32> = Vec::new();
                each_canonical_kmer(bytes, k, |code| {
                    if let Some(h) = sampled(code, cutoff) {
                        sampled_n += 1;
                        let c = *counts.get(&h).unwrap_or(&0);
                        if c >= floor {
                            kept.push(c);
                        }
                    }
                });
                let n_kept = kept.len();
                let mean = if n_kept == 0 {
                    0.0
                } else {
                    kept.iter().map(|&c| c as f64).sum::<f64>() / n_kept as f64
                };
                let med = median(&mut kept);
                (med, mean, gc_fraction(bytes), bytes.len(), sampled_n, n_kept)
            })
            .collect()
    }
}

#[pymethods]
impl DepthIndex {
    /// `k` is the k-mer length, at most 31. `scale` of 1 counts everything,
    /// 100 keeps about one k-mer in a hundred.
    #[new]
    #[pyo3(signature = (k=21, scale=100))]
    pub fn new(k: usize, scale: u64) -> PyResult<Self> {
        if k == 0 || k > 31 {
            return Err(pyo3::exceptions::PyValueError::new_err(
                "k must be between 1 and 31, because a k-mer is packed into a 64 bit integer",
            ));
        }
        if scale == 0 {
            return Err(pyo3::exceptions::PyValueError::new_err(
                "scale must be 1 or more",
            ));
        }
        Ok(DepthIndex {
            k,
            scale,
            cutoff: if scale == 1 {
                u64::MAX
            } else {
                u64::MAX / scale
            },
            counts: FxHashMap::default(),
        })
    }

    #[getter]
    pub fn k(&self) -> usize {
        self.k
    }

    #[getter]
    pub fn scale(&self) -> u64 {
        self.scale
    }

    /// How many distinct sampled k-mers are in the table.
    #[getter]
    pub fn n_kmers(&self) -> usize {
        self.counts.len()
    }

    /// Add a batch of reads to the counts.
    ///
    /// Each worker counts its own share into a private table and the tables are
    /// then folded together, which avoids every thread fighting over one lock.
    pub fn add_reads(&mut self, py: Python<'_>, reads: Vec<String>) {
        py.allow_threads(move || self.count_into(&reads));
    }

    /// Describe each read: median k-mer count, mean k-mer count, GC share,
    /// length, and how many sampled k-mers went into those numbers.
    ///
    /// The median is the one that matters. A mean is pulled around by any short
    /// stretch of repeat inside an otherwise ordinary read, and the whole point
    /// is to describe where the read came from rather than its worst few
    /// hundred bases.
    #[pyo3(signature = (reads, min_count=1))]
    pub fn profile_reads(
        &self,
        py: Python<'_>,
        reads: Vec<String>,
        min_count: u32,
    ) -> Vec<(f64, f64, f64, usize, usize, usize)> {
        py.allow_threads(|| self.profile_into(&reads, min_count))
    }

    /// How many distinct k-mers sit at each count, from 1 up to `max_depth`.
    ///
    /// Anything above the top is piled into the last slot so nothing is lost.
    /// This is the curve the compartment fitting is run on.
    #[pyo3(signature = (max_depth=10000))]
    pub fn depth_histogram(&self, max_depth: usize) -> Vec<u64> {
        let mut hist = vec![0u64; max_depth + 1];
        for &c in self.counts.values() {
            let slot = (c as usize).min(max_depth);
            hist[slot] += 1;
        }
        hist
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn median_handles_both_lengths_and_the_empty_case() {
        assert_eq!(median(&mut [3, 1, 2]), 2.0);
        assert_eq!(median(&mut [4, 1, 3, 2]), 2.5);
        assert_eq!(median(&mut []), 0.0);
    }

    #[test]
    fn sampling_is_a_property_of_the_kmer_not_the_read() {
        // The same k-mer must always land on the same side of the cutoff,
        // otherwise counts from different reads could not be compared.
        let cutoff = u64::MAX / 100;
        for code in 0u64..5000 {
            assert_eq!(sampled(code, cutoff).is_some(), sampled(code, cutoff).is_some());
        }
    }

    #[test]
    fn a_depth_floor_rescues_the_median_when_most_kmers_are_errors() {
        // A short true sequence seen many times, then a long stretch of
        // sequence seen only once. That is what a read looks like to a k-mer
        // counter when most of its k-mers carry a sequencing error: a minority
        // that the rest of the read set agrees with, and a majority that is
        // unique to this read.
        let mut idx = DepthIndex::new(11, 1).unwrap();
        let truth: String = {
            let mut x = 12345u64;
            (0..80)
                .map(|_| {
                    x = crate::kmer::splitmix64(x);
                    [b'A', b'C', b'G', b'T'][(x >> 33) as usize % 4] as char
                })
                .collect()
        };
        let junk: String = {
            let mut x = 99999u64;
            (0..400)
                .map(|_| {
                    x = crate::kmer::splitmix64(x);
                    [b'A', b'C', b'G', b'T'][(x >> 33) as usize % 4] as char
                })
                .collect()
        };

        let copies: Vec<String> = (0..50).map(|_| truth.clone()).collect();
        idx.count_into(&copies);

        let noisy = format!("{truth}{junk}");
        idx.count_into(&[noisy.clone()]);

        let no_floor = idx.profile_into(&[noisy.clone()], 1);
        let with_floor = idx.profile_into(&[noisy], 5);

        // Most k-mers are unique, so the plain median sits down at the junk.
        assert!(
            no_floor[0].0 <= 2.0,
            "without a floor the median should land on the junk, got {}",
            no_floor[0].0
        );
        // With a floor it reports the depth of the sequence that is really shared.
        assert!(
            with_floor[0].0 > 40.0,
            "with a floor the median should report the real depth, got {}",
            with_floor[0].0
        );
        assert!(with_floor[0].5 > 0, "some k-mers should survive the floor");
    }

    #[test]
    fn the_sampled_share_is_close_to_one_over_scale() {
        let scale = 50u64;
        let cutoff = u64::MAX / scale;
        let kept = (0u64..200_000).filter(|&c| sampled(c, cutoff).is_some()).count();
        let expected = 200_000.0 / scale as f64;
        // Allow a generous band. This is a sanity check on the hash, not a
        // statistical test.
        assert!(
            (kept as f64) > expected * 0.7 && (kept as f64) < expected * 1.3,
            "kept {kept}, expected about {expected}"
        );
    }
}
