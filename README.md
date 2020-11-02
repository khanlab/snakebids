
# Snakebids
[![Tests](https://github.com/snakebids/bids_path/actions/workflows/test.yml/badge.svg?branch=main)](https://github.com/snakebids/bids_path/actions/workflows/test.yml?query=branch%3Amain)
[![Version](https://img.shields.io/github/v/tag/snakebids/bids_path?label=version)](https://pypi.org/project/bids_path/)
[![Python versions](https://img.shields.io/pypi/pyversions/bids_path)](https://pypi.org/project/bids_path/)
[![DOI](https://zenodo.org/badge/309495236.svg)](https://zenodo.org/badge/latestdoi/309495236)
 [![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

This repository implements the Snakebids builder function for Bids Paths. You can install this library as a stand-alone package via pip (see below), however, its contents are repackaged in the main [snakebids repository](https://github.com/akhanf/snakebids).

Snakebids is a Python package that extends [Snakemake](https://snakemake.github.io), enabling users to create reproducible, scalable pipelines for processing neuroimaging data in the [BIDS format](https://bids.neuroimaging.io). Snakebids workflows expose a CLI that conforms to the [BIDS App](https://bids-apps.neuroimaging.io) guidelines.

## Features
The `bids()` function is a flexible, unopinionated tool to build valid bids paths in a clean, pythonic syntax. It tries to support all the official bids tags, but allows any arbitrary custom tags as well. Entities are automatically ordered according to the [bids spec](https://bids-specification.readthedocs.io/en/stable/index.html).

## Installation
Snakebids can be installed using pip:

```bash
pip install bids_path
```

## Usage
To use, just import the `bids()` function:

```py
from bids_path import bids
```

Note that `bids()` is just one piece of the broader [`snakebids`](https://github.com/akhanf/snakebids). If `snakebids` is installed, `bids()` can be imported directly from the `snakebids` namespace:

```py
from snakebids import bids
```

For detailed instructions and examples, please refer to the [**documentation**](https://snakebids.readthedocs.io/en/stable/index.html).

## Contributing
Snakebids is an open-source project, and contributions are welcome! If you have any bug reports, feature requests, or improvements, please submit them to the [**issues page**](https://github.com/snakebids/bids_path).

To contribute, first clone the Github repository. `bids_path` dependencies are managed with Poetry (version 1.2 or higher). Please refer to the [poetry website](https://python-poetry.org/docs/master/#installation) for installation instructions.

_Note: `bids_path` makes use of Poetry's dynamic versioning. To see a version number on locally installed Snakebids versions, you will have to also install `poetry-dynamic-versioning` plugin to your poetry installation (`poetry self add "poetry-dynamic-versioning\[plugin\]"). This is **not required** for contribution._

Following installation of Poetry, the development can be set up by running the following commands:

```bash
poetry install
poetry run poe setup
```

`bids_path` uses [poethepoet](https://github.com/nat-n/poethepoet) as a task runner. You can see what commands are available by running:

```bash
poetry run poe
```

Tests are done with `pytest` and can be run via:

```bash
poetry run poe test
```

Additionally, `bids_path` uses pre-commit hooks (installed via the `poe setup` command above) to lint and format code (we use [black](https://github.com/psf/black), [isort](https://github.com/PyCQA/isort) and [ruff](https://beta.ruff.rs/docs/)). By default, these hooks are run on every commit. Please be sure they all pass before making a PR.

## License
`bids_path` is distributed under the MIT License.
