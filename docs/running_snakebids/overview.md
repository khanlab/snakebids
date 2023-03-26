Overview
========

Once you've specified a snakebids app with a config file and one or more workflow files, you're ready to invoke your snakebids app with the standard BIDS app CLI.

Snakebids apps generated with the cookiecutter template will have a simple executable called `run.py`, which exposes the BIDS app CLI, with any additional options configured in the `snakebids.yml` config file. Installing the project with pip will also add an executable with the project's name to the path. Any Snakemake arguments should be added to the end of the invocation.

While Snakebids apps use the standard BIDS app CLI (i.e. `{app_name} {input} {output} {analysis_level}`), it is possible to override the input location for each input type defined in the configuration file. By passing the path override argument `--path_{input_type} {path}`, where `{input_type}` is the name of an input type defined in the configuration file, and `{path}` is the path to a directly containing files appropriate for that input type. If a path override argument is provided for every input type, an `{input}` argument must still be provided to the BIDS app CLI, but it does not need to be an existing directory; a string like `-` will work.

Note that if any rules in the Snakebids workflow use Singularity containers, special precautions must be taken to ensure the input dataset is bound to the Singularity environment. Either:

1. Inputs are copied into a working subdirectory of the output directory before any processing that requires a Singularity container is performed, or:
2. The `SINGULARITY_BINDPATH` environment variable binds the location of the input dataset.

Indexing of large datasets can be a time-consuming process. Snakebids, through `PyBIDS` has the ability to create or leverage an existing database, requiring indexing of datasets to be only performed when user chooses to do so (usually if the dataset has changed)! Note, this feature is **opt-in**, meaning it is not used unless the associated config variables are used. To opt-in:

1. Uncomment the lines in `snakebids.yml` containing `pybids_db_dir` and `pybids_db_reset`.
1. The variables can be updated directly in this file or through the CLI by using `-pybidsdb-dir {dir}` to specify the database path and `--reset-db` to indicate that the database should be updated. _Note: CLI arguments take precendence if both CLI and config variables are set._

Workflow mode
=============

Snakebids apps use a BIDS app CLI, giving great flexibility when switching datasets. However, when developing a Snakebids app or when running the app repeatedly on the same dataset, it can be more convenient to directly call the Snakemake CLI. Snakebids facilitates this using workflow mode.

Workflow mode activates when the Snakebids app itself is used as the output folder. Snakebids will save your config file in the config folder and put any outputs in the results folder. After the first Snakebids call, the Snakemake CLI can be called directly using the generated config file.

As an example, suppose we have a BIDS formatted dataset:

    dataset
    ├── code
    ├── dataset_description.json
    ├── derivatives
    ├── sub-001
    ├── sub-002
    ├── sub-003
    ├── sub-004
    ├── sub-005
    ├── sub-006
    ├── sub-007
    ├── sub-008
    └── sub-009

We'd like to develop a new processing pipeline called SuperCorrect. We'll start by making a new directory in ``derivatives`` using the Snakemake CookieCutter template:

    dataset
    ├── code
    ├── dataset_description.json
    ├── derivatives
    │   └── super_correct
    │       ├── config
    │       │   └── snakebids.yml
    │       ├── pipeline_description.json
    │       ├── run.py
    │       └── workflow
    │           └── Snakefile
    ├── sub-001
    ├── sub-002
    ├── sub-003
    ├── sub-004
    ├── sub-005
    ├── sub-006
    ├── sub-007
    ├── sub-008
    └── sub-009

`cd` to `super_correct`:

    super_correct
    ├── config
    │   └── snakebids.yml
    ├── pipeline_description.json
    ├── run.py
    └── workflow
        └── Snakefile

We then develop our pipeline, writing our `Snakefile`, rules, etc. When we're ready to start testing, we start by calling our `run.py` function:

    run.py ../../ . participant [snakemake args]

We'll see a message telling us the app is running in snakemake mode and, if our workflow doesn't have any bugs, the app will run! `config/snakebids.yml` will be updated to include all the information we passed into the CLI. Output files will be in the `results` folder:

    super_correct
    ├── config
    │   └── snakebids.yml
    ├── pipeline_description.json
    ├── results
    │   └── super_correct
    │       ├── dataset_description.json
    │       ├── sub-001
    │       ├── sub-002
    │       ├── sub-003
    │       ├── sub-004
    │       ├── sub-005
    │       ├── sub-006
    │       ├── sub-007
    │       ├── sub-008
    │       └── sub-009
    ├── run.py
    └── workflow
        └── Snakefile


From now on, instead of calling `run.py`, we can just the Snakemake CLI directly. It will use the same inputs and outputs saved into our config by Snakebids:

    snakemake [args]

You can still use the Snakebids CLI on other datasets. However, if you plan on modifying any files, including config, to make the Snakebids app suitable for the new dataset, it's recommended to use git to clone the app into the `derivatives` folder of the new dataset. Alternatively, you can call ``run.py`` with the `--workflow-mode` flag:

    run.py /path/to/newdata /path/to/newdata/derivatives/super_correct participant --workflow-mode [snakemake args]

This will make a copy of the Snakebids app at the new output directory, excluding the `results` folder and some configuration folders/files (`.snakemake/`, `.snakebids`). It will, again, make a new config file, and put new results in the `output/results` folder.
