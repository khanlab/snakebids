# Snakebids

[![Tests](https://github.com/khanlab/snakebids/actions/workflows/test.yml/badge.svg?branch=main)](https://github.com/khanlab/snakebids/actions/workflows/test.yml?query=branch%3Amain)
[![codecov](https://codecov.io/gh/khanlab/snakebids/graph/badge.svg?token=Q15K5HX7W9)](https://codecov.io/gh/khanlab/snakebids)
[![Documentation Status](https://readthedocs.org/projects/snakebids/badge/?version=stable)](https://snakebids.readthedocs.io/en/stable/?badge=stable)
[![Version](https://img.shields.io/github/v/tag/khanlab/snakebids?label=version)](https://pypi.org/project/snakebids/)
[![Python versions](https://img.shields.io/pypi/pyversions/snakebids)](https://pypi.org/project/snakebids/)
[![DOI](https://zenodo.org/badge/309495236.svg)](https://zenodo.org/badge/latestdoi/309495236)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

Snakebids is a Python package that extends [Snakemake](https://snakemake.github.io), enabling users to create reproducible, scalable pipelines for processing neuroimaging data in the [BIDS format](https://bids.neuroimaging.io). Snakebids workflows expose a CLI that conforms to the [BIDS App](https://bids-apps.neuroimaging.io) guidelines.

## Features

Snakebids includes all of the features of Snakemake, including flexible configuration, parallel execution, and Docker/Singularity support, plus:

- **Built-in support for BIDS datasets**: Seamless workflow functionality with a wide range of BIDS datasets, accommodating various levels of complexity.
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

To contribute, first clone the Github repository.

```bash
git clone https://github.com/khanlab/snakebids
cd snakebids
```

Snakebids dependencies are managed with uv. This can be installed with their standalone installer:

```bash
curl -LsSf https://astral.sh/uv/install.sh | sh
```

> [!TIP]
> Please refer to the [uv website](https://docs.astral.sh/uv/getting-started/installation/) for detailed installation instructions and usage information.

Before coding, run the following command to setup our pre-commit hooks:

```bash
uv run poe setup
```

These check the code against [ruff](https://docs.astral.sh/ruff/) and [pyright](https://github.com/microsoft/pyright). Be sure they both pass before making a PR.


> [!TIP]
> Snakebids uses [poethepoet](https://github.com/nat-n/poethepoet) as a task runner. If this tool is installed globally, it will automatically detect the snakebids environment when directly called. So, instead of `uv run poe [COMMAND]`, you can call `poe [COMMAND]`. It can be installed with:
> ```bash
> uv tool poethepoet
> ```

To check code quality, use:

```bash
uv run poe quality
```

Tests are done with `pytest` and can be run via:

```bash
uv run poe test
```

## License

Snakebids is distributed under the MIT License.

## Acknowledgements

Snakebids extends the Snakemake workflow management system and follows the guidelines outlined by the BIDS specification.

## Relevant papers

- Mölder F, Jablonski KP, Letcher B et al. Sustainable data analysis with Snakemake [version 2; peer review: 2 approved]. F1000Research. 2021. doi: [10.12688/f1000research.29032.2](https://doi.org/10.12688/f1000research.29032.2)
