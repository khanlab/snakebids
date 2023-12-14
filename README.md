# Snakebids

[![Tests](https://github.com/khanlab/snakebids/actions/workflows/test.yml/badge.svg?branch=main)](https://github.com/khanlab/snakebids/actions/workflows/test.yml?query=branch%3Amain)
[![Documentation Status](https://readthedocs.org/projects/snakebids/badge/?version=stable)](https://snakebids.readthedocs.io/en/stable/?badge=stable)
[![Version](https://img.shields.io/github/v/tag/khanlab/snakebids?label=version)](https://pypi.org/project/snakebids/)
[![Python versions](https://img.shields.io/pypi/pyversions/snakebids)](https://pypi.org/project/snakebids/)
[![DOI](https://zenodo.org/badge/309495236.svg)](https://zenodo.org/badge/latestdoi/309495236)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

Snakebids is a Python package that extends [Snakemake](https://snakemake.github.io), enabling users to create reproducible, scalable pipelines for processing neuroimaging data in the [BIDS format](https://bids.neuroimaging.io). Snakebids workflows expose a CLI that conforms to the [BIDS App](https://bids-apps.neuroimaging.io) guidelines.

## Features

Snakebids includes all of the features of Snakemake, including flexible configuration, parallel execution, and Docker/Singularity support, plus:

- **Built-in support for BIDS datasets**: Seamless workflow functionality with a wide range of BIDS datasets, accomodating various levels of complexity.
- **BIDS App Creation**: Provide command-line invocations of your workflow following BIDS App guidelines, ensuring reproducibility and enhancing accessibility of your workflow.
- **BIDS Path Construction**: Easy, flexible construction of valid BIDS paths following BIDS guiding principles, promoting data organization and sharing.
- **Plugin System**: Extend the functionality of Snakebids by creating and using plugins to meet your workflow's needs.
- **Pybids Querying**: Leverages [Pybids](https://bids-standard.github.io/pybids/) to efficiently retrieve specific data required.

## Installation

Snakebids can be installed using pip:

```bash
pip install snakebids
```

## Usage

To create and run a Snakebids workflow, you need to:

1. **Create a Snakefile**: Define the steps of your workflow, including input / output files, processing rules, and dependencies
1. **Create a configuration file**: Customize workflow behaviour using a YAML configuration file. Specify input / output directories and custom workflow parameters.
1. **Run the pipeline**: Execute the Snakebids pipeline by invoking the BIDS App CLI or via Snakemake executable.

For detailed instructions and examples, please refer to the [**documentation**](https://snakebids.readthedocs.io/en/stable/index.html).

## Contributing

Snakebids is an open-source project, and contributions are welcome! If you have any bug reports, feature requests, or improvements, please submit them to the [**issues page**](https://github.com/khanlab/snakebids).

To contribute, first clone the Github repository. Snakebids dependencies are managed with Poetry (version 1.2 or higher). Please refer to the [poetry website](https://python-poetry.org/docs/master/#installation) for installation instructions.

_Note: Snakebids makes use of Poetry's dynamic versioning. To see a version number on locally installed Snakebids versions, you will have to also install `poetry-dynamic-versioning` plugin to your poetry installation (`poetry self add "poetry-dynamic-versioning\[plugin\]"). This is **not required** for contribution._

Following installation of Poetry, the development can be set up by running the following commands:

```bash
poetry install
poetry run poe setup
```

Snakebids uses [poethepoet](https://github.com/nat-n/poethepoet) as a task runner. You can see what commands are available by running:

```bash
poetry run poe
```

Tests are done with `pytest` and can be run via:

```bash
poetry run poe test
```

Additionally, Snakebids uses pre-commit hooks (installed via the `poe setup` command above) to lint and format code (we use [black](https://github.com/psf/black), [isort](https://github.com/PyCQA/isort) and [ruff](https://beta.ruff.rs/docs/). By default, these hooks are run on every commit. Please be sure they all pass before making a PR.

## License

Snakebids is distributed under the MIT License.

## Acknowledgements

Snakebids extends the Snakemake workflow management system and follows the guidelines outlined by the BIDS specification.

## Relevant papers

- Mölder F, Jablonski KP, Letcher B et al. Sustainable data analysis with Snakemake [version 2; peer review: 2 approved]. F1000Research. 2021. doi: [10.12688/f1000research.29032.2](https://doi.org/10.12688/f1000research.29032.2)
