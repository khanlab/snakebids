
snakebids
=========
.. image:: https://readthedocs.org/projects/snakebids/badge/?version=latest
  :target: https://snakebids.readthedocs.io/en/latest/?badge=latest
  :alt: Documentation Status

Snakemake + BIDS

This package allows you to build BIDS Apps using Snakemake. It offers:


* Flexible data grabbing with PyBIDS, configurable solely by config file entries
* Helper function for creating BIDS paths inside Snakemake workflows/rules
* Command-line invocation of snakemake workflows with BIDS App compliance
* Configurable argument parsing specified using the Snakemake workflow config
* Execution either as command-line BIDS apps or via snakemake executable

Contributing
============

Clone the git repository. Snakebids dependencies are managed with Poetry, which you'll need installed on your machine. You can find instructions on the `poetry website <https://python-poetry.org/docs/master/#installation>`_. Then, setup the development environment with the following commands::

  poetry install
  poetry run poe setup

Snakebids uses `poethepoet <https://github.com/nat-n/poethepoet>`_ as a task runner. You can see what commands are available by running::

    poetry run poe

If you wish, you can also run ``poe [[command]]`` directly by installing ``poethepoet`` on your system. Follow the install instructions at the link above.

Tests are done with ``pytest`` and can be run via::

  poetry run pytest

Snakebids uses pre-commit hooks (installed via the ``poe setup`` command above) to lint and format code (we use `black <https://github.com/psf/black>`_, `isort <https://github.com/PyCQA/isort>`_, `pylint <https://pylint.org/>`_ and `flake8 <https://flake8.pycqa.org/en/latest/>`_). By default, these hooks are run on every commit. Please be sure they all pass before making a PR.
