Overview
========

Once you've specified a snakebids app with a config file and one or more workflow files, you're ready to invoke your snakebids app with the standard BIDS app CLI.

Snakebids apps generated with the cookiecutter template will have a simple executable called ``run.py``, which exposes the BIDS app CLI, with any additional options configured in the ``snakebids.yml`` config file. Installing the project with pip will also add an executable with the project's name to the path. Any Snakemake arguments should be added to the end of the invocation.

Note that if any rules in the snakebids workflow use Singularity containers, special precautions must be taken to ensure the input dataset is bound to the Singularity environment. Either:

1. Inputs are copied into a working subdirectory of the output directory before any processing that requires a Singularity container is performed, or:
2. The ``SINGULARITY_BINDPATH`` environment variable binds the location of the input dataset.
