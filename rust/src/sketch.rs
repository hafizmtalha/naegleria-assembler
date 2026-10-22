//! Picking anchor k-mers, and using them to find which reads overlap.
//!
//! Comparing every read against every other read is far too slow, so instead we
//! reduce each read to a short list of anchor k-mers and only compare reads that
//! share some. The usual name for these anchors is minimizers: take a sliding
//! window of k-mers and keep the smallest one in each window. Two reads that
//! overlap will choose the same anchors in the shared stretch, because the choice
//! only depends on the sequence.
//!
//! ## Two changes for an AT rich genome
//!
//! The textbook version orders k-mers alphabetically and keeps the smallest. In a
//! genome that is about two thirds A and T that is a bad idea. Alphabetical order
//! puts A first, so runs of A win nearly every window they appear in. Those runs
//! are everywhere, they say almost nothing about position, and lining two reads up
//! on one is close to meaningless. We order by hash instead, so base composition
//! stops deciding which anchors get picked.
//!
//! On top of that, any k-mer whose own base composition is too lopsided is dropped
//! before the window ever sees it. A run of A, a stretch of AT repeat and similar
//! low complexity sequence all fail that check. What is left are anchors that
//! actually pin down a place in the genome.

use crate::kmer::{base_code, kmer_entropy, splitmix64};
use pyo3::prelude::*;
use rayon::prelude::*;
use rustc_hash::{FxHashMap, FxHashSet};

/// Every canonical k-mer of a sequence as (hash, position, entropy).
fn hashed_kmers(seq: &[u8], k: usize) -> Vec<(u64, u32, f64, bool)> {
    let mut out = Vec::new();
    if k == 0 || k > 31 || seq.len() < k {
        return out;
    }
    let mask: u64 = (1u64 << (2 * k)) - 1;
    let shift = 2 * (k as u64 - 1);
    let mut fwd: u64 = 0;
    let mut rev: u64 = 0;
    let mut run = 0usize;

    for (i, &b) in seq.iter().enumerate() {
        match base_code(b) {
            Some(c) => {
                fwd = ((fwd << 2) | c) & mask;
                rev = (rev >> 2) | ((3 - c) << shift);
                run += 1;
                if run >= k {
                    let forward = fwd <= rev;
                    let canon = if forward { fwd } else { rev };
                    let pos = (i + 1 - k) as u32;
                    out.push((splitmix64(canon), pos, kmer_entropy(canon, k), forward));
                }
            }
            None => {
                run = 0;
                fwd = 0;
                rev = 0;
            }
        }
    }
    out
}

/// Anchor k-mers for one sequence, as (hash, position) sorted by position.
///
/// `w` is the window width in k-mers. `min_entropy` is in bits, from 0 to 2, and
/// a k-mer scoring below it is never chosen.
pub fn minimizers(seq: &[u8], k: usize, w: usize, min_entropy: f64) -> Vec<(u64, u32, bool)> {
    let all = hashed_kmers(seq, k);
    let usable: Vec<(u64, u32, bool)> = all
        .into_iter()
        .filter(|&(_, _, e, _)| e >= min_entropy)
        .map(|(h, p, _, f)| (h, p, f))
        .collect();

    if usable.is_empty() || w == 0 {
        return usable;
    }
    if usable.len() <= w {
        // Too short for a full window. Keep the single best anchor.
        let best = usable.iter().min_by_key(|&&(h, _, _)| h).copied();
        return best.into_iter().collect();
    }

    // Monotonic deque holding indices of candidates that could still win a
    // window. Anything with a larger hash than a later arrival can be dropped,
    // because the later one survives longer and beats it.
    let mut out: Vec<(u64, u32, bool)> = Vec::new();
    let mut deque: std::collections::VecDeque<usize> = std::collections::VecDeque::new();
    let mut last_picked: Option<u32> = None;

    for i in 0..usable.len() {
        while let Some(&back) = deque.back() {
            if usable[back].0 >= usable[i].0 {
                deque.pop_back();
            } else {
                break;
            }
        }
        deque.push_back(i);

        // Drop anything that has fallen out of the left edge of the window.
        while let Some(&front) = deque.front() {
            if i >= w && front + w <= i {
                deque.pop_front();
            } else {
                break;
            }
        }

        if i + 1 >= w {
            let &front = deque.front().unwrap();
            let pick = usable[front];
            // The same k-mer wins several windows in a row. Record it once.
            if last_picked != Some(pick.1) {
                out.push(pick);
                last_picked = Some(pick.1);
            }
        }
    }
    out
}

/// Anchor k-mers for one sequence. Exposed so the Python side and the tests can
/// look at what the picker does.
#[pyfunction]
#[pyo3(signature = (seq, k=15, w=10, min_entropy=1.0))]
pub fn sketch(seq: &str, k: usize, w: usize, min_entropy: f64) -> Vec<(u64, u32, bool)> {
    minimizers(seq.as_bytes(), k, w, min_entropy)
}

/// Holds the anchors of a set of reads and reports which pairs share enough of
/// them to be worth aligning properly.
#[pyclass]
pub struct OverlapIndex {
    k: usize,
    w: usize,
    min_entropy: f64,
    sketches: Vec<Vec<(u64, u32, bool)>>,
}

impl OverlapIndex {
    /// Sketch a batch of reads, with no Python involved.
    pub fn sketch_into(&mut self, reads: &[String]) {
        let (k, w, e) = (self.k, self.w, self.min_entropy);
        let mut new: Vec<Vec<(u64, u32, bool)>> = reads
            .par_iter()
            .map(|r| minimizers(r.as_bytes(), k, w, e))
            .collect();
        self.sketches.append(&mut new);
    }

    /// The pair search itself, with no Python involved.
    pub fn overlaps_into(
        &self,
        min_shared: usize,
        max_occurrence: usize,
    ) -> Vec<(usize, usize, usize)> {
        // Which reads carry each anchor.
        let mut by_anchor: FxHashMap<u64, Vec<u32>> = FxHashMap::default();
        for (read_id, sk) in self.sketches.iter().enumerate() {
            let mut seen: FxHashSet<u64> = FxHashSet::default();
            for &(h, _, _) in sk {
                if seen.insert(h) {
                    by_anchor.entry(h).or_default().push(read_id as u32);
                }
            }
        }
        by_anchor.retain(|_, reads| reads.len() <= max_occurrence && reads.len() > 1);

        // For each read, count how many anchors it shares with each partner that
        // comes later in the list. Doing only the later ones means each pair is
        // counted once and no deduplication is needed afterwards.
        (0..self.sketches.len())
            .into_par_iter()
            .flat_map_iter(|a| {
                let mut shared: FxHashMap<u32, usize> = FxHashMap::default();
                let mut seen: FxHashSet<u64> = FxHashSet::default();
                for &(h, _, _) in &self.sketches[a] {
                    if !seen.insert(h) {
                        continue;
                    }
                    if let Some(reads) = by_anchor.get(&h) {
                        for &b in reads {
                            if b as usize > a {
                                *shared.entry(b).or_insert(0) += 1;
                            }
                        }
                    }
                }
                shared
                    .into_iter()
                    .filter(move |&(_, n)| n >= min_shared)
                    .map(move |(b, n)| (a, b as usize, n))
                    .collect::<Vec<_>>()
            })
            .collect()
    }
}

/// Where two reads line up, worked out from the anchors they share.
#[derive(Clone, Copy, Debug)]
pub struct OverlapHit {
    pub a: usize,
    pub b: usize,
    pub same_strand: bool,
    pub shift: i64,
    pub anchors: usize,
    pub a_start: u32,
    pub a_end: u32,
    pub b_start: u32,
    pub b_end: u32,
}

impl OverlapIndex {
    /// Turn a candidate pair into an actual placement, or decide there isn't one.
    ///
    /// Two reads that genuinely overlap share anchors whose positions move in
    /// step. Subtract one position from the other for every shared anchor and
    /// the real ones all land on the same value, give or take indels. Anchors
    /// that match by chance scatter. So we bin those differences and keep the
    /// busiest bin. That bin is the shift between the two reads, and the anchors
    /// in it say how far the agreement runs.
    ///
    /// When the reads come from opposite strands the positions run against each
    /// other, so their sum stays constant instead of their difference. Anchors
    /// are canonical, meaning a k-mer and its reverse complement hash the same,
    /// which is what lets a flipped pair match at all. The stored strand flag
    /// says which of the two cases a given anchor match is.
    fn place_pair(
        &self,
        a: usize,
        b: usize,
        bin_width: i64,
        min_anchors: usize,
    ) -> Option<OverlapHit> {
        let sa = &self.sketches[a];
        let sb = &self.sketches[b];
        if sa.is_empty() || sb.is_empty() {
            return None;
        }

        let mut by_hash: FxHashMap<u64, Vec<(u32, bool)>> = FxHashMap::default();
        for &(h, p, f) in sb {
            by_hash.entry(h).or_default().push((p, f));
        }

        // (diagonal value, position in a, position in b, same strand)
        let mut matches: Vec<(i64, u32, u32, bool)> = Vec::new();
        for &(h, pa, fa) in sa {
            if let Some(hits) = by_hash.get(&h) {
                for &(pb, fb) in hits {
                    let same = fa == fb;
                    let diag = if same {
                        pa as i64 - pb as i64
                    } else {
                        pa as i64 + pb as i64
                    };
                    matches.push((diag, pa, pb, same));
                }
            }
        }
        if matches.len() < min_anchors {
            return None;
        }

        // Bin the diagonals. Indels shift a diagonal a little, so a bin a few
        // tens of bases wide keeps a real chain together without merging
        // chains that are genuinely different.
        let width = bin_width.max(1);
        let mut bins: FxHashMap<(i64, bool), Vec<usize>> = FxHashMap::default();
        for (i, &(diag, _, _, same)) in matches.iter().enumerate() {
            bins.entry((diag.div_euclid(width), same)).or_default().push(i);
        }

        // A chain can straddle a bin edge, so each bin is scored together with
        // its right hand neighbour.
        let mut best: Option<((i64, bool), Vec<usize>)> = None;
        let keys: Vec<(i64, bool)> = bins.keys().copied().collect();
        for key in keys {
            let mut members = bins.get(&key).cloned().unwrap_or_default();
            if let Some(next) = bins.get(&(key.0 + 1, key.1)) {
                members.extend_from_slice(next);
            }
            if best.as_ref().map_or(true, |(_, m)| members.len() > m.len()) {
                best = Some((key, members));
            }
        }

        let (key, members) = best?;
        if members.len() < min_anchors {
            return None;
        }

        let same_strand = key.1;
        let mut a_start = u32::MAX;
        let mut a_end = 0u32;
        let mut b_start = u32::MAX;
        let mut b_end = 0u32;
        let mut diag_sum: i64 = 0;
        for &i in &members {
            let (diag, pa, pb, _) = matches[i];
            diag_sum += diag;
            a_start = a_start.min(pa);
            a_end = a_end.max(pa);
            b_start = b_start.min(pb);
            b_end = b_end.max(pb);
        }

        Some(OverlapHit {
            a,
            b,
            same_strand,
            shift: diag_sum / members.len() as i64,
            anchors: members.len(),
            a_start,
            a_end,
            b_start,
            b_end,
        })
    }

    /// Place every candidate pair. This is the plain Rust entry point.
    pub fn place_all(
        &self,
        min_shared: usize,
        max_occurrence: usize,
        bin_width: i64,
        min_anchors: usize,
    ) -> Vec<OverlapHit> {
        let candidates = self.overlaps_into(min_shared, max_occurrence);
        candidates
            .into_par_iter()
            .filter_map(|(a, b, _)| self.place_pair(a, b, bin_width, min_anchors))
            .collect()
    }
}

#[pymethods]
impl OverlapIndex {
    #[new]
    #[pyo3(signature = (k=15, w=10, min_entropy=1.0))]
    pub fn new(k: usize, w: usize, min_entropy: f64) -> PyResult<Self> {
        if k == 0 || k > 31 {
            return Err(pyo3::exceptions::PyValueError::new_err(
                "k must be between 1 and 31",
            ));
        }
        Ok(OverlapIndex {
            k,
            w,
            min_entropy,
            sketches: Vec::new(),
        })
    }

    /// Sketch and store a batch of reads. Reads keep the order they arrive in,
    /// and that order is what the pair indices refer to later.
    pub fn add_reads(&mut self, py: Python<'_>, reads: Vec<String>) {
        py.allow_threads(move || self.sketch_into(&reads));
    }

    #[getter]
    pub fn n_reads(&self) -> usize {
        self.sketches.len()
    }

    /// How many anchors each read ended up with.
    pub fn sketch_sizes(&self) -> Vec<usize> {
        self.sketches.iter().map(|s| s.len()).collect()
    }

    /// Candidate overlapping pairs, as (read a, read b, anchors shared) with
    /// a < b.
    ///
    /// `min_shared` is how many anchors two reads must have in common before the
    /// pair is worth reporting.
    ///
    /// `max_occurrence` throws away any anchor that turns up in more than that
    /// many reads. Those come from repeats, and keeping them would both bury the
    /// real signal and make the work blow up, since an anchor in n reads costs n
    /// squared comparisons on its own.
    /// Place candidate pairs and give back where each pair lines up, as
    /// (read a, read b, same strand, shift, anchors, a start, a end, b start, b end).
    #[pyo3(signature = (min_shared=3, max_occurrence=500, bin_width=64, min_anchors=4))]
    pub fn place_overlaps(
        &self,
        py: Python<'_>,
        min_shared: usize,
        max_occurrence: usize,
        bin_width: i64,
        min_anchors: usize,
    ) -> Vec<(usize, usize, bool, i64, usize, u32, u32, u32, u32)> {
        py.allow_threads(|| {
            self.place_all(min_shared, max_occurrence, bin_width, min_anchors)
                .into_iter()
                .map(|h| {
                    (
                        h.a,
                        h.b,
                        h.same_strand,
                        h.shift,
                        h.anchors,
                        h.a_start,
                        h.a_end,
                        h.b_start,
                        h.b_end,
                    )
                })
                .collect()
        })
    }

    #[pyo3(signature = (min_shared=3, max_occurrence=500))]
    pub fn find_overlaps(
        &self,
        py: Python<'_>,
        min_shared: usize,
        max_occurrence: usize,
    ) -> Vec<(usize, usize, usize)> {
        py.allow_threads(|| self.overlaps_into(min_shared, max_occurrence))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn anchors_sit_inside_the_sequence() {
        let seq = b"ACGTTGCAAGGTCCATGACTGCATTTACGGCATCAGTACCGTAAC";
        let m = minimizers(seq, 11, 5, 1.0);
        assert!(!m.is_empty());
        for &(_, pos, _) in &m {
            assert!(pos as usize + 11 <= seq.len());
        }
    }

    #[test]
    fn anchors_come_out_in_position_order() {
        let seq = b"ACGTTGCAAGGTCCATGACTGCATTTACGGCATCAGTACCGTAACGGATCAGT";
        let m = minimizers(seq, 11, 5, 1.0);
        for pair in m.windows(2) {
            assert!(pair[0].1 < pair[1].1);
        }
    }

    #[test]
    fn a_poly_a_run_contributes_no_anchors() {
        // The entropy filter is the whole point of this test. Without it these
        // k-mers would win every window they touch.
        let m = minimizers(&[b'A'; 200], 15, 10, 1.0);
        assert!(m.is_empty(), "poly-A should never be chosen as an anchor");
    }

    #[test]
    fn dropping_the_entropy_filter_lets_poly_a_back_in() {
        let m = minimizers(&[b'A'; 200], 15, 10, 0.0);
        assert!(!m.is_empty());
    }

    fn random_dna(n: usize, seed: u64) -> Vec<u8> {
        let mut x = seed | 1;
        let bases = [b'A', b'C', b'G', b'T'];
        (0..n)
            .map(|_| {
                x = crate::kmer::splitmix64(x);
                bases[(x >> 33) as usize % 4]
            })
            .collect()
    }

    #[test]
    fn placement_recovers_a_known_shift() {
        let whole = random_dna(4000, 7);
        // b starts 1200 bases into a, so a's coordinates are 1200 ahead of b's.
        let a = &whole[0..3000];
        let b = &whole[1200..4000];

        let mut idx = OverlapIndex::new(15, 10, 1.0).unwrap();
        idx.sketch_into(&[
            String::from_utf8(a.to_vec()).unwrap(),
            String::from_utf8(b.to_vec()).unwrap(),
        ]);

        let hit = idx.place_pair(0, 1, 64, 4).expect("the pair should be placed");
        assert!(hit.same_strand);
        assert!(
            (hit.shift - 1200).abs() <= 64,
            "shift was {} but should be about 1200",
            hit.shift
        );
        assert!(hit.anchors >= 4);
    }

    #[test]
    fn placement_finds_a_flipped_pair() {
        let whole = random_dna(4000, 11);
        let a = whole[0..3000].to_vec();
        let mut b: Vec<u8> = whole[1200..4000]
            .iter()
            .rev()
            .map(|&c| match c {
                b'A' => b'T',
                b'T' => b'A',
                b'C' => b'G',
                _ => b'C',
            })
            .collect();
        b.truncate(2800);

        let mut idx = OverlapIndex::new(15, 10, 1.0).unwrap();
        idx.sketch_into(&[
            String::from_utf8(a).unwrap(),
            String::from_utf8(b).unwrap(),
        ]);

        let hit = idx.place_pair(0, 1, 64, 4).expect("the flipped pair should be placed");
        assert!(!hit.same_strand, "this pair is reverse complemented");
        assert!(hit.anchors >= 4);
    }

    #[test]
    fn unrelated_sequences_are_not_placed() {
        let a = random_dna(3000, 3);
        let b = random_dna(3000, 99);
        let mut idx = OverlapIndex::new(15, 10, 1.0).unwrap();
        idx.sketch_into(&[
            String::from_utf8(a).unwrap(),
            String::from_utf8(b).unwrap(),
        ]);
        assert!(idx.place_pair(0, 1, 64, 4).is_none());
    }

    #[test]
    fn overlapping_sequences_share_anchors_and_unrelated_ones_do_not() {
        // One stretch of DNA read twice with a shared middle.
        let whole: Vec<u8> = b"GCATTACGGATCCAGTTACGGCATCAGTACCGTAACGGATCAGTTGCACCATGGACTTAGCCATGGACTTAGC"
            .to_vec();
        let a = &whole[..55];
        let b = &whole[20..];
        let other = b"TTTTGGGGCCCCAAAATTTTGGGGCCCCAAAATTTTGGGGCCCCAAAATTTTGGGG";

        let mut idx = OverlapIndex::new(11, 5, 1.0).unwrap();
        idx.sketches.push(minimizers(a, 11, 5, 1.0));
        idx.sketches.push(minimizers(b, 11, 5, 1.0));
        idx.sketches.push(minimizers(other, 11, 5, 1.0));

        let hits = idx.overlaps_into(1, 500);
        assert!(
            hits.iter().any(|&(x, y, _)| x == 0 && y == 1),
            "the two halves of one sequence should be found"
        );
        assert!(
            !hits.iter().any(|&(_, y, _)| y == 2),
            "an unrelated sequence should not be paired with anything"
        );
    }
}
