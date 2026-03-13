# Optional Rust Acceleration

`snakebids` ships with an optional compiled extension
(`snakebids._rust._core`) that accelerates `SnakemakeFormatter.parse()`.
The extension is a private implementation detail—no public API changes
are required to use it.

## How it works

When the extension is present, `SnakemakeFormatter.parse()` delegates
its inner loop to a Rust implementation built with
[PyO3](https://pyo3.rs/). The Python fallback is always available and is
used automatically when the compiled extension is absent.

## Building the extension locally

You will need [Rust](https://rustup.rs/) and
[maturin](https://www.maturin.rs/) installed.

```bash
# install maturin once
pip install maturin

# build and place the extension inside the source tree
maturin develop --release
```

After this, `snakebids._rust._core` is importable and all
`SnakemakeFormatter` calls benefit from the faster implementation
transparently.

## Verifying the extension is active

```python
from snakebids.utils.snakemake_templates import _HAS_RUST_PARSE
print(_HAS_RUST_PARSE)  # True when the extension is loaded
```

## CI / packaging notes

* Pure-Python installs (`pip install snakebids`) **do not** require
  Rust—the extension is optional.
* Pre-built wheels distributed on PyPI will include the compiled
  extension for supported platforms.
* The `[tool.maturin]` section in `pyproject.toml` configures the
  module placement (`snakebids._rust._core`) and Python source root.
