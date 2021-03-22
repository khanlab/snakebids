Overview
========

Once you've specified a snakebids app with a config file and one or more workflow files, you're ready to invoke your snakebids app with the standard BIDS app CLI.

Snakebids apps generated with the cookiecutter template will have a simple executable called ``run.py``, which exposes the BIDS app CLI, with any additional options configured in the ``snakebids.yml`` config file. Installing the project with pip will also add an executable with the project's name to the path. Any Snakemake arguments should be added to the end of the invocation.

While Snakebids apps use the standard BIDS app CLI (i.e. ``{app_name} {input} {output} {analysis_level}``), it is possible to override the input location for each input type defined in the configuration file. By passing the path override argument ``--path_{input_type} {path}``, where ``{input_type}`` is the name of an input type defined in the configuration file, and ``{path}`` is the path to a directly containing files appropriate for that input type. If a path override argument is provided for every input type, an ``{input}`` argument must still be provided to the BIDS app CLI, but it does not need to be an existing directory; a string like ``-`` will work.

Note that if any rules in the snakebids workflow use Singularity containers, special precautions must be taken to ensure the input dataset is bound to the Singularity environment. Either:

1. Inputs are copied into a working subdirectory of the output directory before any processing that requires a Singularity container is performed, or:
2. The ``SINGULARITY_BINDPATH`` environment variable binds the location of the input dataset.
