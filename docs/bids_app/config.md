{#bids-app-config}
# Configuration

Snakebids is configured with a YAML (or JSON) file that extends the standard [snakemake config file](https://snakemake.readthedocs.io/en/stable/snakefiles/configuration.html#standard-configuration) with variables that snakebids uses to parse an input BIDS dataset and expose the snakebids workflow to the command line.

## Config Variables

### `pybids_inputs`

A dictionary that describes each type of input you want to grab from an input BIDS dataset. Snakebids will parse your dataset with {func}`generate_inputs() <snakebids.generate_inputs>`, converting each input type into a {class}`BidsComponent <snakebids.BidsComponent>`. The value of each item should be a dictionary with keys `filters` and `wildcards`.

#### Filters

The value of `filters` should be a dictionary where each key corresponds to a BIDS entity, and the value specifies which values of that entity should be grabbed. The dictionary for each input is sent to the [PyBIDS' `get()` function ](#bids.layout.BIDSLayout). `filters` can be set according to a few different formats:

* [`string`](#str): specifies an exact value for the entity. In the following example:
  ```yaml
  pybids_inputs:
    bold:
      filters:
        suffix: 'bold'
        extension: '.nii.gz'
        datatype: 'func'
  ```

  the bold component would match any paths under the `func/` datatype folder, with the suffix `bold` and the extension `.nii.gz`.

  ```
  sub-xxx/.../func/sub-xxx_ses-xxx_..._bold.nii.gz
  ```

* [`boolean`](#bool): constrains presence or absence of the entity without restricting its value. `False` requires that the entity be **absent**, while `True` requires the entity to be **present**, regardless of value.
  ```yaml
  pybids_inputs:
    derivs:
      filters:
        datatype: 'func'
        desc: True
        acquisition: False
  ```
  The above example selects all paths in the `func/` datatype folder that have a `_desc-` entity but do not have the `_acq-` entity.

* [`list`](#list): Specify multiple string or boolean filters. Any path matching any one of the filters will be selected. Using `False` as one of the filters allows the entity to optionally be absent in addition to matching one of the string filters. Using `True` along with text is redundant, as `True` will cause any value to be selected. Using `True` with `False` is equivalent to not providing the filter at all.

  These filters:

  ```yaml
  pybids_inputs:
    derivs:
      filters:
        acquisition:
          - False
          - MPRAGE
          - MP2RAGE
  ```

  would select all of the following paths:

  ```
  sub-001/ses-1/anat/sub-001_ses-001_acq-MPRAGE_run-1_T1w.nii.gz
  sub-001/ses-1/anat/sub-001_ses-001_acq-MP2RAGE_run-1_T1w.nii.gz
  sub-001/ses-1/anat/sub-001_ses-001_run-1_T1w.nii.gz
  ```


* To use regex for filtering, use an additional subkey set either to [`match`](#re.match) or [`search`](#re.search), depending on which regex method you wish to use. This key may be set to any one of the above items (`str`, `bool`, or `list`). Only one such key may be used.

  These filters:

  ```yaml
  pybids_inputs:
    derivs:
      filters:
        suffix:
          search: '[Tt]1'
        acquisition:
          match: MP2?RAGE
  ```

  would select all of the following paths:

  ```
  sub-001/ses-1/anat/sub-001_ses-001_acq-MPRAGE_run-1_T1.nii.gz
  sub-001/ses-1/anat/sub-001_ses-001_acq-MP2RAGE_run-1_t1w.nii.gz
  sub-001/ses-1/anat/sub-001_ses-001_acq-MPRAGE_run-1_qT1w.nii.gz
  ```

````{note}
`match` and `search` are both _filtering methods_. In addition to these, `get` is also a valid filtering method and may be used as the subkey for a filter. However, this is equivalent to directly providing the desired filter without a subkey:

```yaml
pybids_inputs:
  derivs:
    filters:
      suffix:
        get: T1w

# is the same as
pybids_inputs:
  derivs:
    filters:
      suffix: T1w
```

In other words, `get` is the default filtering method.
````

#### Wildcards

The value of `wildcards` should be a list of BIDS entities. Snakebids collects the values of any entities specified and saves them in the {attr}`entities <snakebids.BidsComponent.entities>` and {attr}`~snakebids.BidsComponent.zip_lists` entries of the corresponding {class}`BidsComponent <snakebids.BidsComponent>`. In other words, these are the entities to be preserved in output paths derived from the input being described. Placing an entity in `wildcards` does not require the entity be present. If an entity is not found, it will be left out of {attr}`entities <snakebids.BidsComponent.entities>`. To require the presence of an entity, place it under `filters` set to `true`.

In the following (YAML-formatted) example, the `bold` input type is specified. BIDS files with the datatype `func`, suffix `bold`, and extension `.nii.gz` will be grabbed, and the `subject`, `session`, `acquisition`, `task`, and `run` entities of those files will be left as wildcards. The `task` entity must be present, but there must not be any `desc`.

```yaml
pybids_inputs:
  bold:
    filters:
      suffix: 'bold'
      extension: '.nii.gz'
      datatype: 'func'
      task: true
      desc: false
    wildcards:
      - subject
      - session
      - acquisition
      - task
      - run
```

### `pybidsdb_dir`

PyBIDS allows for the use of a cached layout to be used in order to reduce the time required to index a BIDS dataset. A path (if provided) to save the *pybids* [layout](#bids.layout.BIDSLayout). If `None` or `''` is provided, the layout is not saved or used. The path provided must be absolute, otherwise the database will not be used.

### `pybidsdb_reset`

A boolean determining whether the existing layout should be be updated. Default behaviour does not update the existing database if one is used.

### `analysis_levels`

A list of analysis levels in the BIDS app. Typically, this will include participant and/or group. Note that the default (YAML) configuration file expects this mapping to be identified with the anchor ``analysis_levels`` to be aliased by ``parse_args``.

### `targets_by_analysis_level`

A mapping from the name of each ``analysis_level`` to the list of rules or files to be run for that analysis level.

### `parse_args`

A dictionary of command-line parameters to make available as part of the BIDS app. Each item of the mapping is passed to [argparse's `add_argument` function](#argparse.ArgumentParser.add_argument). A number of default entries are present in a new snakebids project's config file that structure the BIDS app's CLI, but additional command-line arguments can be added as necessary.

As in [`ArgumentParser.add_argument()`](#argparse.ArgumentParser.add_argument), `type` may be used to convert the argument to the specified type. It may be set to any type that can be serialized into yaml, for instance, `str`, `int`, `float`, and `boolean`.

```yaml
parse_args:
  --a-string:
    help: args are string by default
  --a-path:
    help: |
      A path pointing to data needed for the pipeline. These are still converted
      into strings, but are first resolved into absolute paths (see below)
    type: Path
  --another-path:
    help: This type annotation does the same thing as above
    type: pathlib.Path
  --a-number:
    help: A number important for the analysis
    type: float
```

When CLI parameters are used to collect paths, `type` should be set to [`Path`](#pathlib.Path) (or [`pathlib.Path`](#pathlib.Path)). These arguments will still be serialized as strings (since yaml doesn't have a path type), but snakebids will automatically resolve all arguments into absolute paths. This is important to prevent issues with snakebids and relative paths.


### `debug`

A boolean that determines whether debug statements are printed during parsing. Should be disabled (False) if you're generating DAG visualization with snakemake.


### `derivatives`

A boolean (or path(s) to derivatives datasets) that determines whether snakebids will search in the derivatives subdirectory of the input dataset.
