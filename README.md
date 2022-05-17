
snakebids
=========
[![Tests](https://github.com/akhanf/snakebids/actions/workflows/test.yml/badge.svg?branch=main)](https://github.com/akhanf/snakebids/actions/workflows/test.yml?query=branch%3Amain)
[![Documentation Status](https://readthedocs.org/projects/snakebids/badge/?version=stable)](https://snakebids.readthedocs.io/en/stable/?badge=stable)
[![Version](https://img.shields.io/github/v/tag/akhanf/snakebids?label=version)](https://pypi.org/project/snakebids/)
[![Python versions](https://img.shields.io/pypi/pyversions/snakebids)](https://pypi.org/project/snakebids/)
[![DOI](https://zenodo.org/badge/309495236.svg)](https://zenodo.org/badge/latestdoi/309495236)

Snakemake + BIDS
This package allows you to build BIDS Apps using Snakemake. It offers:


* Flexible data grabbing with PyBIDS, configurable solely by config file entries
* Helper function for creating BIDS paths inside Snakemake workflows/rules
* Command-line invocation of snakemake workflows with BIDS App compliance
* Configurable argument parsing specified using the Snakemake workflow config
* Execution either as command-line BIDS apps or via snakemake executable

Contributing
============

Clone the git repository. Snakebids dependencies are managed with Poetry, which you'll need installed on your machine. You can find instructions on the [poetry website](https://python-poetry.org/docs/master/#installation). Then, setup the development environment with the following commands:

```bash
poetry install
poetry run poe setup
```

Snakebids uses [poethepoet](https://github.com/nat-n/poethepoet) as a task runner. You can see what commands are available by running:

```bash
poetry run poe
```

If you wish, you can also run `poe [command]` directly by installing `poethepoet` on your system. Follow the install instructions at the link above.

Tests are done with `pytest` and can be run via:

```bash
poetry run pytest
```

Snakebids uses pre-commit hooks (installed via the `poe setup` command above) to lint and format code (we use [black](https://github.com/psf/black), [isort](https://github.com/PyCQA/isort), [pylint](https://pylint.org/) and [flake8](https://flake8.pycqa.org/en/latest/)). By default, these hooks are run on every commit. Please be sure they all pass before making a PR.
