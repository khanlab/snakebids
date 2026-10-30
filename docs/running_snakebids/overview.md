(running-snakebids)=
# Running Snakebids

Once you've specified a snakebids app with a config file and one or more workflow files, you're ready to invoke your snakebids app with the standard BIDS app CLI.

Snakebids apps generated with the cookiecutter template will have a simple executable called `run.py`, which exposes the BIDS app CLI, with any additional options configured in the `snakebids.yml` config file. Installing the project with pip will also add an executable with the project's name to the path. Any Snakemake arguments should be added to the end of the invocation.

While Snakebids apps use the standard BIDS app CLI (i.e. `{app_name} {input} {output} {analysis_level}`), it is possible to override the input location for each input type defined in the configuration file. By passing the path override argument `--path_{input_type} {path}`, where `{input_type}` is the name of an input type defined in the configuration file, and `{path}` is the path to a directly containing files appropriate for that input type. If a path override argument is provided for every input type, an `{input}` argument must still be provided to the BIDS app CLI, but it does not need to be an existing directory; a string like `-` will work.

Note that if any rules in the Snakebids workflow use Singularity containers, special precautions must be taken to ensure the input dataset is bound to the Singularity environment. Either:

1. Inputs are copied into a working subdirectory of the output directory before any processing that requires a Singularity container is performed, or:
2. The `SINGULARITY_BINDPATH` environment variable binds the location of the input dataset.

Indexing of large datasets can be a time-consuming process. Leveraging the functionality of `PyBIDS`, Snakebids offers a convenient solution by allowing you to create or utilize an existing database. With this approach, the indexing of datasets is only performed when explicitly requested, typically following changes to the dataset. To create or use an existing database, you can invoke the following CLI arguments:

1. `--pybidsdb-dir {dir}`: specify the path to the database directory
1. `--pybidsdb-reset`: indicate that an existing database should be updated

The boilerplate app starts with the validator plugin enabled - without it, validation is not performed. By default, this feature uses the command-line (node.js) version of the [validator](https://www.npmjs.com/package/bids-validator). If this is not found to be installed on the system, the `pybids` version of validation will be performed instead. To opt-out of validation, invoke the `--skip-bids-validation` flag. Details related to using and creating plugins can be found on the [plugins](/bids_app/plugins) page.

(selecting-participants)=
## Selecting participants

By default, every subject found in the dataset is processed. Three arguments narrow this down:

- `--participant-label LABEL [LABEL ...]`: only process the listed subjects
- `--exclude-participant-label LABEL [LABEL ...]`: process every subject except those listed
- `--participant-filter QUERY`: only process the subjects whose row in the dataset's `participants.tsv` matches `QUERY`

Labels are given without the `sub-` prefix (e.g. `001`, not `sub-001`).

### Filtering on `participants.tsv`

`--participant-filter` selects subjects using the columns of the [`participants.tsv`](https://bids-specification.readthedocs.io/en/stable/modality-agnostic-files.html#participants-file) file at the root of the dataset. Given the following file:

```text
participant_id	age	sex	group
sub-001	34	F	control
sub-002	71	M	patient
sub-003	n/a	F	patient
```

the following call processes only `sub-002`:

```console
$ my_app bids/ out/ participant --participant-filter "age > 60 and group == 'patient'"
```

The query is evaluated with {meth}`pandas.DataFrame.query`, so it supports comparisons (`==`, `!=`, `<`, `>=`, ...), `and`/`or`/`not`, membership tests (`group in ['patient', 'prodromal']`), and string methods (`group.str.startswith('pat')`). Some practical notes:

- Wrap the whole query in double quotes and use single quotes for text values, so the shell passes it to the app unchanged.
- The subject label (without `sub-`) is available as `subject`, e.g. `"subject in ['001', '005'] or age < 30"`.
- Column names containing spaces or characters such as `-` must be surrounded by backticks: ``"`handedness-score` > 0"``.
- `n/a` values are read as missing, and never satisfy a comparison: in the example above, `sub-003` matches neither `age > 60` nor `age <= 60`. Use `age.isna()` to select subjects with missing values.
- Column types are inferred from their contents. A numeric column containing any non-numeric value other than `n/a` is read as text, and numeric comparisons on it will fail.

`--participant-filter` can be combined with the other two arguments: subjects must match the query *and* be listed in `--participant-label` (if given), and must not be listed in `--exclude-participant-label`.

Subjects that have data in the dataset but no row in `participants.tsv` are excluded when a filter is given, and a warning lists them. A warning is also given if the query matches no subjects. If the dataset has no `participants.tsv`, or the query cannot be evaluated (for instance, because it names a column that does not exist), the app stops with an error listing the available columns.

```{note}
`--participant-filter` is only available in apps whose Snakefile passes `participant_filter=config.get("participant_filter")` to {func}`~snakebids.generate_inputs`. Apps created from the snakebids template include this line; older apps need to add it (see [](/bids_app/workflow)).
```

### Limitations

- **Components with a `subject` filter in their configuration are not filtered.** If an input component in the app's config already filters on `subject`, the participant selection is ignored for that component. This also applies to `--participant-label` and `--exclude-participant-label`.
- **No matches with custom paths.** If a query matches no subjects and a component uses a custom path (`--path-{component}`), the app stops with an error (`Boolean filters in items with custom paths are not supported`) instead of returning an empty component. The warning that no participants matched is shown just before this error.
- **Only the dataset root is used.** `participants.tsv` is read from the root of the input dataset. The `participants.tsv` of derivative datasets, and the column descriptions in `participants.json`, are not used.
- **Participant-level metadata only.** Values that change between sessions, such as those recorded in `sessions.tsv` (e.g. age at each scan), cannot be used for filtering.
- **The query is code.** {meth}`pandas.DataFrame.query` evaluates the query as a Python expression. Treat `--participant-filter` like any other command-line argument, and don't pass it text from untrusted sources.

## Workflow mode

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
    │   └── super_correct
    │       ├── config
    │       │   └── snakebids.yml
    │       ├── pipeline_description.json
    │       ├── run.py
    │       └── workflow
    │           └── Snakefile
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
    │   └── snakebids.yml
    ├── pipeline_description.json
    ├── run.py
    └── workflow
        └── Snakefile

We then develop our pipeline, writing our `Snakefile`, rules, etc. When we're ready to start testing, we start by calling our `run.py` function:

    run.py ../../ . participant [snakemake args]

We'll see a message telling us the app is running in snakemake mode and, if our workflow doesn't have any bugs, the app will run! `config/snakebids.yml` will be updated to include all the information we passed into the CLI. Output files will be in the `results` folder:

    super_correct
    ├── config
    │   └── snakebids.yml
    ├── pipeline_description.json
    ├── results
    │   └── super_correct
    │       ├── dataset_description.json
    │       ├── sub-001
    │       ├── sub-002
    │       ├── sub-003
    │       ├── sub-004
    │       ├── sub-005
    │       ├── sub-006
    │       ├── sub-007
    │       ├── sub-008
    │       └── sub-009
    ├── run.py
    └── workflow
        └── Snakefile


From now on, instead of calling `run.py`, we can just the Snakemake CLI directly. It will use the same inputs and outputs saved into our config by Snakebids:

    snakemake [args]

You can still use the Snakebids CLI on other datasets. However, if you plan on modifying any files, including config, to make the Snakebids app suitable for the new dataset, it's recommended to use git to clone the app into the `derivatives` folder of the new dataset. Alternatively, you can call ``run.py`` with the `--workflow-mode` flag:

    run.py /path/to/newdata /path/to/newdata/derivatives/super_correct participant --workflow-mode [snakemake args]

This will make a copy of the Snakebids app at the new output directory, excluding the `results` folder and some configuration folders/files (`.snakemake/`, `.snakebids`). It will, again, make a new config file, and put new results in the `output/results` folder.
