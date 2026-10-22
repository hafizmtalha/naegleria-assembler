//! Turning DNA text into numbers, and a few things that follow from that.
//!
//! Every k-mer up to length 31 fits in a u64 at two bits per base, so the whole
//! library works on integers rather than strings. A k-mer and its reverse
//! complement describe the same piece of DNA, so we always keep whichever of the
//! two is the smaller number. That is what "canonical" means below.

/// A, C, G and T map to 0 to 3. Anything else, including N, gives None and
/// breaks the current run of k-mers.
#[inline]
pub fn base_code(b: u8) -> Option<u64> {
    match b {
        b'A' | b'a' => Some(0),
        b'C' | b'c' => Some(1),
        b'G' | b'g' => Some(2),
        b'T' | b't' => Some(3),
        _ => None,
    }
}

/// splitmix64. A cheap way to spread bit patterns out evenly.
///
/// We need this because the raw two bit encoding is not random at all. Runs of A
/// encode to small numbers and runs of T to large ones, so picking k-mers by
/// their raw value would pick by base composition. Hashing first removes that.
#[inline]
pub fn splitmix64(x: u64) -> u64 {
    let mut z = x.wrapping_add(0x9E37_79B9_7F4A_7C15);
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    z ^ (z >> 31)
}

/// Walk a sequence and hand every canonical k-mer to the closure.
///
/// The k-mer is carried forward one base at a time in both directions at once,
/// so each new base costs a couple of shifts rather than a rescan.
#[inline]
pub fn each_canonical_kmer<F: FnMut(u64)>(seq: &[u8], k: usize, mut f: F) {
    if k == 0 || k > 31 || seq.len() < k {
        return;
    }
    let mask: u64 = (1u64 << (2 * k)) - 1;
    let shift = 2 * (k as u64 - 1);
    let mut fwd: u64 = 0;
    let mut rev: u64 = 0;
    let mut run = 0usize;

    for &b in seq {
        match base_code(b) {
            Some(c) => {
                fwd = ((fwd << 2) | c) & mask;
                rev = (rev >> 2) | ((3 - c) << shift);
                run += 1;
                if run >= k {
                    f(if fwd < rev { fwd } else { rev });
                }
            }
            None => {
                // An N means we cannot trust the window that spans it.
                run = 0;
                fwd = 0;
                rev = 0;
            }
        }
    }
}

/// Shannon entropy of the base composition inside one encoded k-mer, in bits.
///
/// A k-mer of a single repeated base scores 0. An even mix of all four scores 2.
/// This is the cheapest useful way to spot the low complexity stretches that an
/// AT rich genome is full of.
#[inline]
pub fn kmer_entropy(code: u64, k: usize) -> f64 {
    let mut counts = [0u32; 4];
    let mut c = code;
    for _ in 0..k {
        counts[(c & 3) as usize] += 1;
        c >>= 2;
    }
    let n = k as f64;
    let mut h = 0.0f64;
    for &cnt in counts.iter() {
        if cnt > 0 {
            let p = cnt as f64 / n;
            h -= p * p.log2();
        }
    }
    h
}

/// Collapse runs of the same base down to one.
///
/// Long read platforms lose count inside homopolymers, so "AAAAA" and "AAAA" are
/// the same stretch of DNA read twice with a different answer. Comparing reads on
/// the collapsed view makes that disagreement disappear. The original sequence is
/// kept elsewhere, because the collapsed form is only for finding overlaps, never
/// for the final answer.
pub fn homopolymer_compress(seq: &[u8]) -> (Vec<u8>, Vec<u32>) {
    let mut out = Vec::with_capacity(seq.len());
    let mut runs: Vec<u32> = Vec::with_capacity(seq.len());
    let mut last: Option<u8> = None;
    for &b in seq {
        let up = b.to_ascii_uppercase();
        if Some(up) == last {
            if let Some(r) = runs.last_mut() {
                *r += 1;
            }
        } else {
            out.push(up);
            runs.push(1);
            last = Some(up);
        }
    }
    (out, runs)
}

/// Share of G and C among the bases we can read. N and friends are ignored.
pub fn gc_fraction(seq: &[u8]) -> f64 {
    let mut gc = 0usize;
    let mut known = 0usize;
    for &b in seq {
        match b {
            b'G' | b'g' | b'C' | b'c' => {
                gc += 1;
                known += 1;
            }
            b'A' | b'a' | b'T' | b't' => known += 1,
            _ => {}
        }
    }
    if known == 0 {
        0.0
    } else {
        gc as f64 / known as f64
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn collect(seq: &str, k: usize) -> Vec<u64> {
        let mut v = Vec::new();
        each_canonical_kmer(seq.as_bytes(), k, |x| v.push(x));
        v
    }

    #[test]
    fn a_sequence_and_its_reverse_complement_give_the_same_kmers() {
        let fwd = collect("ACGTTGCAAGGT", 5);
        let mut rev = collect("ACCTTGCAACGT", 5);
        rev.reverse();
        assert_eq!(fwd, rev);
    }

    #[test]
    fn an_n_breaks_the_run() {
        // 10 bases, k = 4, so 7 k-mers if the N were not there.
        assert_eq!(collect("ACGTACGTAC", 4).len(), 7);
        // The N kills every window that touches it.
        assert_eq!(collect("ACGTNCGTAC", 4).len(), 3);
    }

    #[test]
    fn short_sequences_give_nothing() {
        assert!(collect("ACG", 5).is_empty());
    }

    #[test]
    fn entropy_is_zero_for_one_repeated_base_and_two_for_an_even_mix() {
        let mut poly = Vec::new();
        each_canonical_kmer(b"AAAAAAAAAA", 8, |x| poly.push(x));
        assert!(kmer_entropy(poly[0], 8) < 1e-9);

        let mut mixed = Vec::new();
        each_canonical_kmer(b"ACGTACGT", 8, |x| mixed.push(x));
        assert!((kmer_entropy(mixed[0], 8) - 2.0).abs() < 1e-9);
    }

    #[test]
    fn homopolymers_collapse_and_keep_their_lengths() {
        let (seq, runs) = homopolymer_compress(b"AAACCGTTTT");
        assert_eq!(seq, b"ACGT");
        assert_eq!(runs, vec![3, 2, 1, 4]);
    }

    #[test]
    fn gc_ignores_unknown_bases() {
        assert!((gc_fraction(b"GCAT") - 0.5).abs() < 1e-9);
        assert!((gc_fraction(b"GCATNNNN") - 0.5).abs() < 1e-9);
        assert_eq!(gc_fraction(b"NNNN"), 0.0);
    }
}
