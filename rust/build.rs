//! Tells the test binary where to find libpython at run time.
//!
//! The real build is a Python extension module. The host interpreter is already
//! in memory when it loads, so the library is there and nothing needs to say
//! where it lives.
//!
//! `cargo test` is a different story. It builds an ordinary executable, and that
//! executable has to find libpython by itself. Without help it does not, and the
//! tests die before the first one runs:
//!
//! ```text
//! dyld: Library not loaded: @rpath/libpython3.12.dylib
//! ```
//!
//! The compile is fine. Every test passes once the binary can start. It just
//! cannot start. So this asks the interpreter that is building us where its
//! library directory is, and writes that path into the binary as an rpath. Then
//! `cargo test --no-default-features` works on its own, with nothing to set up
//! first, which is the whole point.
//!
//! Only used for the test build. The extension module build skips it, because
//! the `extension-module` feature is on by default and that feature means pyo3
//! does not link libpython at all.

use std::process::Command;

fn main() {
    println!("cargo:rerun-if-changed=build.rs");

    // With the default features on we are building the extension module, which
    // does not link libpython. Nothing to do.
    if std::env::var("CARGO_FEATURE_EXTENSION_MODULE").is_ok() {
        return;
    }

    // Ask whichever python is on PATH where it keeps its libraries. This is the
    // same interpreter maturin would use, so it is the right one.
    let python = std::env::var("PYO3_PYTHON").unwrap_or_else(|_| "python3".to_string());
    let output = Command::new(&python)
        .args(["-c", "import sysconfig; print(sysconfig.get_config_var('LIBDIR') or '')"])
        .output();

    let libdir = match output {
        Ok(out) if out.status.success() => String::from_utf8_lossy(&out.stdout).trim().to_string(),
        // No python, or it would not answer. Say so and carry on rather than
        // failing the build, since the compile itself does not need this.
        _ => {
            println!(
                "cargo:warning=could not ask {python} where libpython lives, so the test \
                 binary may not start. Set DYLD_FALLBACK_LIBRARY_PATH (macOS) or \
                 LD_LIBRARY_PATH (Linux) to the lib folder of your environment."
            );
            return;
        }
    };

    if libdir.is_empty() {
        return;
    }

    println!("cargo:rustc-link-search=native={libdir}");
    println!("cargo:rustc-link-arg=-Wl,-rpath,{libdir}");
}
