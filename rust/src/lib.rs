//! nagcore: the parts of the Naegleria assembler that have to be fast.
//!
//! Only the inner loops live here. Counting k-mers across gigabytes of reads,
//! sketching them and matching the sketches are the steps where Python would
//! spend all day. Everything that decides what any of it means is in Python,
//! where it can be read and changed.

use pyo3::prelude::*;

pub mod depth;
pub mod kmer;
pub mod sketch;

/// Collapse runs of the same base. Gives back the collapsed sequence and the
/// length of each run, so the original can be rebuilt.
#[pyfunction]
fn homopolymer_compress(seq: &str) -> (String, Vec<u32>) {
    let (bytes, runs) = kmer::homopolymer_compress(seq.as_bytes());
    (String::from_utf8_lossy(&bytes).into_owned(), runs)
}

/// Share of G and C among the readable bases.
#[pyfunction]
fn gc_fraction(seq: &str) -> f64 {
    kmer::gc_fraction(seq.as_bytes())
}

#[pymodule]
fn nagcore(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<depth::DepthIndex>()?;
    m.add_class::<sketch::OverlapIndex>()?;
    m.add_function(wrap_pyfunction!(sketch::sketch, m)?)?;
    m.add_function(wrap_pyfunction!(homopolymer_compress, m)?)?;
    m.add_function(wrap_pyfunction!(gc_fraction, m)?)?;
    m.add("__version__", env!("CARGO_PKG_VERSION"))?;
    Ok(())
}
