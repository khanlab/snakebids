"""Private package for optional Rust-accelerated internals.

The compiled extension ``_core`` is built separately with maturin.
If the extension is not available, the pure-Python fallback is used.
"""
