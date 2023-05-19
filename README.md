
# Snakebids
[![Tests](https://github.com/akhanf/snakebids/actions/workflows/test.yml/badge.svg?branch=main)](https://github.com/akhanf/snakebids/actions/workflows/test.yml?query=branch%3Amain)
[![Documentation Status](https://readthedocs.org/projects/snakebids/badge/?version=stable)](https://snakebids.readthedocs.io/en/stable/?badge=stable)
[![Version](https://img.shields.io/github/v/tag/akhanf/snakebids?label=version)](https://pypi.org/project/snakebids/)
[![Python versions](https://img.shields.io/pypi/pyversions/snakebids)](https://pypi.org/project/snakebids/)
[![DOI](https://zenodo.org/badge/309495236.svg)](https://zenodo.org/badge/latestdoi/309495236)
 [![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

Snakebids is a Python package that facilitates the creation of reproducible and scalable data processing pipelines for neuroimaging data in the BIDS format (e.g. BIDS Apps). It leverages the power of Snakemake, a workflow management system, to automate the execution of complex pipelines.

## Features

* **BIDS Compliance**: Helper function to ensure BIDS paths adhere to the BIDS specification, promoting data organization and sharing. Additionally, provide command-line invocations of your workflow with BIDS App compliance.
* **Flexible Configuration**: Easily configure and customize your workflow using YAML configuration files.
* **Plugin System**: Extend the functionality of Snakebids by creating and using plugins.
* **Parallel Execution**: Leverage the power of parallel processing to speed up your workflow.
* **Docker and Singularity Support**: Run your workflow in containerized environments, eliminating the need for installation of required software dependencies and improve reproducibility.

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
Snakebids is an open-source project, and contributions are welcome! If you have any bug reports, feature requests, or improvements, please submit them to the [**issues page**](https://github.com/akhanf/snakebids).

To contribute, first clone the Github repository. Snakebids dependencies are managed with Poetry (version 1.2 or higher). Please refer to the [poetry website](https://python-poetry.org/docs/master/#installation) for installation instructions. Following installation of Poetry, the development can be setup by running the following commands:

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

Additionally, Snakebids uses pre-commit hooks (installed via the `poe setup` command above) to lint and format code (we use [black](https://github.com/psf/black), [isort](https://github.com/PyCQA/isort), [pylint](https://pylint.org/) and [ruff](https://beta.ruff.rs/docs/). By default, these hooks are run on every commit. Please be sure they all pass before making a PR.

## License

Snakebids is distributed under the MIT License.

## Relevant papers
* Mölder F, Jablonski KP, Letcher B et al. Sustainable data analysis with Snakemake [version 2; peer review: 2 approved]. F1000Research. 2021. doi: [10.12688/f1000research.29032.2](https://doi.org/10.12688/f1000research.29032.2){:target="_blank"} 