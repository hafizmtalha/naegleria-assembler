"""nagasm: a genome assembler built around the way Naegleria keeps its DNA.

Most assemblers assume coverage is roughly even across a genome. Naegleria
breaks that assumption badly, because it keeps its ribosomal RNA genes on a small
circle held in thousands of copies, alongside a mitochondrial genome also in many
copies. Those elements are tiny in length and huge in read count.

So this package sorts reads into compartments by depth before assembling
anything, then gives each compartment the treatment it deserves. See
``docs/DESIGN.md`` for the reasoning, and the ``triage`` module for the core of
it.
"""

from __future__ import annotations

__version__ = "0.1.0"

from . import nagcore  # noqa: F401  (the Rust core, built by maturin)

__all__ = [
    "__version__",
    "nagcore",
    "seqio",
    "triage",
]
