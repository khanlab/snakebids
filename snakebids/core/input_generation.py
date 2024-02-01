"""Utilities for converting Snakemake apps to BIDS apps."""
from __future__ import annotations

import json
import logging
import os
import re
import warnings
from collections import defaultdict
from pathlib import Path
from typing import (
    Any,
    Iterable,
    Literal,
    overload,
)

import more_itertools as itx
from bids import BIDSLayout, BIDSLayoutIndexer
from snakemake.script import Snakemake

from snakebids.core._querying import (
    FilterSpecError,
    PostFilter,
    UnifiedFilter,
    get_matching_files,
)
from snakebids.core.datasets import BidsComponent, BidsDataset, BidsDatasetDict
from snakebids.core.filtering import filter_list
from snakebids.exceptions import (
    ConfigError,
    DuplicateComponentError,
    RunError,
)
from snakebids.types import InputConfig, InputsConfig, ZipList
from snakebids.utils.snakemake_io import glob_wildcards
from snakebids.utils.utils import (
    DEPRECATION_FLAG,
    BidsEntity,
    BidsParseError,
    get_first_dir,
)

_logger = logging.getLogger(__name__)


@overload
def generate_inputs(
    bids_dir: Path | str,
    pybids_inputs: InputsConfig,
    pybidsdb_dir: Path | str | None = ...,
    pybidsdb_reset: bool = ...,
    derivatives: bool | Path | str = ...,
    pybids_config: str | None = ...,
    limit_to: Iterable[str] | None = ...,
    participant_label: Iterable[str] | str | None = ...,
    exclude_participant_label: Iterable[str] | str | None = ...,
    use_bids_inputs: Literal[True] | None = ...,
    index_metadata: bool = ...,
    validate: bool = ...,
    pybids_database_dir: Path | str | None = ...,
    pybids_reset_database: bool = ...,
) -> BidsDataset:
    ...


@overload
def generate_inputs(
    bids_dir: Path | str,
    pybids_inputs: InputsConfig,
    pybidsdb_dir: Path | str | None = ...,
    pybidsdb_reset: bool = ...,
    derivatives: bool | Path | str = ...,
    pybids_config: str | None = ...,
    limit_to: Iterable[str] | None = ...,
    participant_label: Iterable[str] | str | None = ...,
    exclude_participant_label: Iterable[str] | str | None = ...,
    use_bids_inputs: Literal[False] = ...,
    index_metadata: bool = ...,
    validate: bool = ...,
    pybids_database_dir: Path | str | None = ...,
    pybids_reset_database: bool = ...,
) -> BidsDatasetDict:
    ...


def generate_inputs(
    bids_dir: Path | str,
    pybids_inputs: InputsConfig,
    pybidsdb_dir: Path | str | None = None,
    pybidsdb_reset: bool | None = None,
    derivatives: bool | Path | str = False,
    pybids_config: str | None = None,
    limit_to: Iterable[str] | None = None,
    participant_label: Iterable[str] | str | None = None,
    exclude_participant_label: Iterable[str] | str | None = None,
    use_bids_inputs: bool | None = None,
    index_metadata: bool = False,
    validate: bool = False,
    pybids_database_dir: Path | str | None = None,
    pybids_reset_database: bool | None = None,
) -> BidsDataset | BidsDatasetDict:
    """Dynamically generate snakemake inputs using pybids_inputs.

    Pybids is used to parse the bids_dir. Custom paths can also be parsed by including
    the custom_paths entry under the pybids_inputs descriptor.

    Parameters
    ----------
    bids_dir
        Path to bids directory

    pybids_inputs
        Configuration for bids inputs, with keys as the names (``str``)

        Nested `dicts` with the following required keys (for complete info, see
        :class:`~snakebids.types.InputConfig`):

        * ``"filters"``: Dictionary of entity: "values" (dict of str -> str or list of
          str). The entity keywords should the bids tags on which to filter. The values
          should be an acceptable str value for that entity, or a list of acceptable str
          values.

        * ``"wildcards"``: List of (str) bids tags to include as wildcards in
          snakemake. At minimum this should usually include
          ``['subject','session']``, plus any other wildcards that you may
          want to make use of in your snakemake workflow, or want to retain
          in the output paths. Any wildcards in this list that are not in the
          filename will just be ignored.

        * ``"custom_path"``: Custom path to be parsed with wildcards wrapped in braces,
          as in ``/path/to/sub-{subject}/{wildcard_1}-{wildcard_2}``. This path will be
          parsed without pybids, allowing the use of non-bids-compliant paths.

    pybidsdb_dir
        Path to database directory. If None is provided, database
        is not used

    pybidsdb_reset
        A boolean that determines whether to reset / overwrite
        existing database.

    derivatives
        Indicates whether pybids should look for derivative datasets under bids_dir.
        These datasets must be properly formatted according to bids specs to be
        recognized. Defaults to False.

    limit_to
        If provided, indicates which input descriptors from pybids_inputs should be
        parsed. For example, if pybids_inputs describes ``"bold"`` and ``"dwi"`` inputs,
        and ``limit_to = ["bold"]``, only the "bold" inputs will be parsed. "dwi" will
        be ignored

    participant_label
        Indicate one or more participants to be included from input parsing. This may
        cause errors if subject filters are also specified in pybids_inputs. It may not
        be specified if exclude_participant_label is specified

    exclude_participant_label
        Indicate one or more participants to be excluded from input parsing. This may
        cause errors if subject filters are also specified in pybids_inputs. It may not
        be specified if participant_label is specified

    use_bids_inputs
        If False, returns the classic :class:`BidsDatasetDict` instead of
        :class`BidsDataset`. Setting to True is deprecated as of v0.8, as this is now
        the default behaviour

    index_metadata
        If True indexes metadata of BIDS directory using pybids, otherwise skips
        indexing.

    validate
        If True performs validation of BIDS directory using pybids, otherwise
        skips validation.

    Returns
    -------
    BidsDataset | BidsDatasetDict
        Object containing organized information about the bids inputs for consumption in
        snakemake. See the documentation of :class:`BidsDataset` for details and
        examples.

    Example
    -------
    As an example, consider the following BIDS dataset::

        example
        ├── README.md
        ├── dataset_description.json
        ├── participant.tsv
        ├── sub-001
        │   ├── ses-01
        │   │   ├── anat
        │   │   │   ├── sub-001_ses-01_run-01_T1w.json
        │   │   │   ├── sub-001_ses-01_run-01_T1w.nii.gz
        │   │   │   ├── sub-001_ses-01_run-02_T1w.json
        │   │   │   └── sub-001_ses-01_run-02_T1w.nii.gz
        │   │   └── func
        │   │       ├── sub-001_ses-01_task-nback_bold.json
        │   │       ├── sub-001_ses-01_task-nback_bold.nii.gz
        │   │       ├── sub-001_ses-01_task-rest_bold.json
        │   │       └── sub-001_ses-01_task-rest_bold.nii.gz
        │   └── ses-02
        │       ├── anat
        │       │   ├── sub-001_ses-02_run-01_T1w.json
        │       │   └── sub-001_ses-02_run-01_T1w.nii.gz
        │       └── func
        │           ├── sub-001_ses-02_task-nback_bold.json
        │           ├── sub-001_ses-02_task-nback_bold.nii.gz
        │           ├── sub-001_ses-02_task-rest_bold.json
        │           └── sub-001_ses-02_task-rest_bold.nii.gz
        └── sub-002
            ├── ses-01
            │   ├── anat
            │   │   ├── sub-002_ses-01_run-01_T1w.json
            │   │   ├── sub-002_ses-01_run-01_T1w.nii.gz
            │   │   ├── sub-002_ses-01_run-02_T1w.json
            │   │   └── sub-002_ses-01_run-02_T1w.nii.gz
            │   └── func
            │       ├── sub-002_ses-01_task-nback_bold.json
            │       ├── sub-002_ses-01_task-nback_bold.nii.gz
            │       ├── sub-002_ses-01_task-rest_bold.json
            │       └── sub-002_ses-01_task-rest_bold.nii.gz
            └── ses-02
                └── anat
                    ├── sub-002_ses-02_run-01_T1w.json
                    ├── sub-002_ses-02_run-01_T1w.nii.gz
                    ├── sub-002_ses-02_run-02_T1w.json
                    └── sub-002_ses-02_run-02_T1w.nii.gz

    With the following ``pybids_inputs`` defined in the config file::

        pybids_inputs:
          bold:
            filters:
              suffix: 'bold'
              extension: '.nii.gz'
              datatype: 'func'
            wildcards:
              - subject
              - session
              - acquisition
              - task
              - run

    Then ``generate_inputs(bids_dir, pybids_input)`` would return the following values::

        BidsDataset({
            "bold": BidsComponent(
                name="bold",
                path="bids/sub-{subject}/ses-{session}/func/sub-{subject}_ses-{session}\
_task-{task}_bold.nii.gz",
                zip_lists={
                    "subject": ["001",   "001",  "001",   "001",  "002",   "002" ],
                    "session": ["01",    "01",   "02",    "02",   "01",    "01"  ],
                    "task":    ["nback", "rest", "nback", "rest", "nback", "rest"],
                },
            ),
            "t1w": BidsComponent(
                name="t1w",
                path="example/sub-{subject}/ses-{session}/anat/sub-{subject}_\
ses-{session}_run-{run}_T1w.nii.gz",
                zip_lists={
                    "subject": ["001", "001", "001", "002", "002", "002", "002"],
                    "session": ["01",  "01",  "02",  "01",  "01",  "02",  "02" ],
                    "run":     ["01",  "02",  "01",  "01",  "02",  "01",  "02" ],
                },
            ),
        })
    """
    postfilters = PostFilter()
    postfilters.add_filter("subject", participant_label, exclude_participant_label)

    pybidsdb_dir, pybidsdb_reset = _normalize_database_args(
        pybidsdb_dir, pybidsdb_reset, pybids_database_dir, pybids_reset_database
    )

    # Generates a BIDSLayout
    layout = (
        _gen_bids_layout(
            bids_dir=bids_dir,
            derivatives=derivatives,
            pybids_config=pybids_config,
            pybidsdb_dir=pybidsdb_dir,
            pybidsdb_reset=pybidsdb_reset,
            index_metadata=index_metadata,
            validate=validate,
        )
        if not _all_custom_paths(pybids_inputs)
        else None
    )

    bids_inputs = _get_components(
        bids_layout=layout,
        inputs_config=pybids_inputs,
        limit_to=limit_to,
        postfilters=postfilters,
    )

    if use_bids_inputs is True:
        _logger.warning(
            "The parameter `use_bids_inputs` in generate_inputs() is now set, by "
            "default, to True. Manually setting it to True is deprecated as of version "
            "0.8. "
        )
    elif use_bids_inputs is None:
        use_bids_inputs = True

    try:
        dataset = BidsDataset.from_iterable(bids_inputs, layout)
    except DuplicateComponentError as err:
        msg = (
            f"Multiple components found with the same name: {err.duplicated_names_str}"
        )
        raise ConfigError(msg) from err

    if use_bids_inputs:
        return dataset
    return dataset.as_dict


def _normalize_database_args(
    pybidsdb_dir: Path | str | None,
    pybidsdb_reset: bool | str | None,
    pybids_database_dir: Path | str | None,
    pybids_reset_database: str | bool | None,
) -> tuple[Path | str | None, bool]:
    """Handle deprecated arguments for pybidsdb."""
    if pybids_database_dir is not None:
        warnings.warn(
            "The parameter `pybids_database_dir` in generate_inputs() is deprecated "
            "and will be removed in the next release. To set the pybids database, use "
            "the `pybidsdb_dir` parameter instead.",
            stacklevel=1,
        )
    if pybids_reset_database is not None:
        warnings.warn(
            "The parameter `pybids_reset_database` in generate_inputs() is deprecated "
            "and will be removed in the next release. To reset the pybids database, "
            "use the `pybidsdb_reset` parameter instead.",
            stacklevel=1,
        )

    pybidsdb_dir = pybidsdb_dir or pybids_database_dir
    pybidsdb_reset = (
        pybidsdb_reset
        if pybidsdb_reset is not None
        else pybids_reset_database
        if pybids_reset_database is not None
        else False
    )

    depr_len = len(DEPRECATION_FLAG)
    if (
        pybidsdb_dir is not None
        and str(pybidsdb_dir).startswith(DEPRECATION_FLAG)
        and str(pybidsdb_dir).endswith(DEPRECATION_FLAG)
    ):
        warnings.warn(
            "The config value `pybids_db_dir` is deprecated and will be removed in a "
            "future release. To access a CLI-specified pybids database directory, use "
            '`config.get("pybidsdb_dir")` instead.',
            stacklevel=1,
        )
        pybidsdb_dir = str(pybidsdb_dir)[depr_len:-depr_len]
    try:
        if (
            isinstance(pybidsdb_reset, str)
            and pybidsdb_reset.startswith(DEPRECATION_FLAG)
            and pybidsdb_reset.endswith(DEPRECATION_FLAG)
        ):
            pybidsdb_reset = bool(int(pybidsdb_reset[depr_len:-depr_len]))
            warnings.warn(
                "The config value `pybids_db_reset` is deprecated and will be removed "
                "in a future release. To access CLI-specified pybids database reset "
                'instructions, use `config.get("pybidsdb_reset")` instead.',
                stacklevel=1,
            )
    except ValueError as err:
        msg = "pybidsdb_reset must be a boolean"
        raise TypeError(msg) from err

    if not isinstance(pybidsdb_reset, bool):
        msg = "pybidsdb_reset must be a boolean"
        raise TypeError(msg)

    return pybidsdb_dir, pybidsdb_reset


def _all_custom_paths(config: InputsConfig):
    """Check that all input components have a custom path."""
    return all(comp.get("custom_path") for comp in config.values())


def _gen_bids_layout(
    *,
    bids_dir: Path | str,
    derivatives: Path | str | bool,
    pybidsdb_dir: Path | str | None,
    pybidsdb_reset: bool,
    pybids_config: Path | str | None = None,
    index_metadata: bool = False,
    validate: bool = False,
) -> BIDSLayout:
    """Create (or reindex) the BIDSLayout.

    Parameters
    ----------
    bids_dir
        Path to bids directory

    derivatives
        A boolean (or path(s) to derivatives datasets) that
        determines whether snakebids will search in the
        derivatives subdirectory of the input dataset.

    pybidsdb_dir
        Path to database directory. If None is provided, database
        is not used

    pybidsdb_reset
        A boolean that determines whether to reset / overwrite
        existing database.

    index_metadata
        A boolen that determines whether to parse and index metadata

    validate
        A boolean that determines whether to validate the bids dataset

    Returns
    -------
    layout : BIDSLayout
        Layout from pybids for accessing the BIDS dataset to grab paths
    """
    # Check for database_dir
    # If blank, assume db not to be used
    if not pybidsdb_dir:
        pybidsdb_dir = None
    # Otherwise check for relative path and update
    elif not Path(pybidsdb_dir).is_absolute():
        pybidsdb_dir = None
        _logger.warning("Absolute path must be provided, database will not be used")

    return BIDSLayout(
        str(bids_dir),
        derivatives=derivatives,
        validate=validate,
        config=pybids_config,
        database_path=pybidsdb_dir,
        reset_database=pybidsdb_reset,
        indexer=BIDSLayoutIndexer(validate=False, index_metadata=index_metadata),
    )


def write_derivative_json(snakemake: Snakemake, **kwargs: dict[str, Any]) -> None:
    """Update sidecar file with provided sources and parameters.

    Intended for usage in snakemake scripts.

    Parameters
    ----------
    snakemake : struct Snakemake
        structure passed to snakemake python scripts, containing input,
        output, params, config ...
        This function requires input.json and output.json to be defined, as
        it will read and write json files
    """
    sidecar = json.loads(Path(snakemake.input.json).read_text())

    sidecar.update(
        {
            "Sources": [snakemake.input],
            "Parameters": snakemake.params,
            **kwargs,
        }
    )

    with open(snakemake.output.json, "w", encoding="utf-8") as outfile:
        json.dump(sidecar, outfile, indent=4)


def _get_components(
    *,
    bids_layout: BIDSLayout | None,
    inputs_config: InputsConfig,
    postfilters: PostFilter,
    limit_to: Iterable[str] | None = None,
):
    """Generate components based on components config and a bids layout.

    Parameters
    ----------
    bids_layout : BIDSLayout
        Layout from pybids for accessing the BIDS dataset to grab paths.

    inputs_config
        Dictionary indexed by modality name, specifying the filters and
        wildcards for each pybids input.

    limit_to
        List of inputs to skip, this used by snakebids to exclude modalities based on
        cmd-line args.

    postfilters
        Filters to all components after delineation.

    Yields
    ------
    BidsComponent:
        One BidsComponent is yielded for each modality described by ``pybids_inputs``.

    Raises
    ------
    ConfigError
        In response to invalid configuration, missing components, or parsing errors.
    """
    for name in limit_to or inputs_config:
        comp = _get_component(
            bids_layout=bids_layout,
            component=inputs_config[name],
            input_name=name,
            postfilters=postfilters,
        )
        if comp is not None:
            yield comp


def _get_component(
    bids_layout: BIDSLayout | None,
    component: InputConfig,
    *,
    input_name: str,
    postfilters: PostFilter,
) -> BidsComponent | None:
    """Create component based on provided config.

    Parameters
    ----------
    bids_layout
        Layout from pybids for accessing the BIDS dataset to grab paths

    component
        Dictionary indexed by modality name, specifying the filters and
        wildcards for each pybids input.

    input_name
        Name of the component.

    postfilters
        Filters to component after delineation

    Raises
    ------
    ConfigError
        In response to invalid configuration, missing components, or parsing errors.
    """
    _logger.debug("Grabbing inputs for %s...", input_name)

    filters = UnifiedFilter(component, postfilters or {})

    if "custom_path" in component:
        path = component["custom_path"]
        zip_lists = _parse_custom_path(path, filters=filters)
        return BidsComponent(name=input_name, path=path, zip_lists=zip_lists)

    if bids_layout is None:
        msg = (
            f"No valid bids dir given, but {input_name} does not have a "
            "custom_path specified."
        )
        raise RunError(msg)

    zip_lists: dict[str, list[str]] = defaultdict(list)
    paths: set[str] = set()
    try:
        matching_files = get_matching_files(bids_layout, filters)
    except FilterSpecError as err:
        raise err.get_config_error(input_name) from err

    for img in matching_files:
        wildcards: list[str] = [
            wildcard
            for wildcard in component.get("wildcards", [])
            if wildcard in img.entities
        ]
        _logger.debug("Wildcards %s found entities for %s", wildcards, img.path)

        try:
            path, parsed_wildcards = _parse_bids_path(img.path, wildcards)
        except BidsParseError as err:
            msg = (
                "Parsing failed:\n"
                f"  Entity: {err.entity.entity}\n"
                f"  Pattern: {err.entity.regex}\n"
                f"  Path: {img.path}\n"
                "\n"
                "Pybids parsed this path using the pattern: "
                f"{bids_layout.entities[err.entity.entity].regex}\n"
                "\n"
                "Snakebids is not currently able to handle this entity. If it is a "
                "custom entity, its `tag-` must be configured to be the same as "
                "its name. Its entry in your pybids config file should look like:\n"
                f'{{\n\t"name": "{err.entity.entity}",\n'
                f'\t"pattern":"{err.entity.entity}-(<value_pattern>)"\n}}\n'
                f"If {err.entity.entity} is an official pybids entity, please "
                "ensure you are using the latest version of snakebids"
            )
            raise ConfigError(msg) from err

        for wildcard_name, value in parsed_wildcards.items():
            zip_lists[wildcard_name].append(value)

        paths.add(path)

    # now, check to see if unique
    if len(paths) == 0:
        _logger.warning(
            "No input files found for snakebids component %s:\n"
            "    filters:\n%s\n"
            "    wildcards:\n%s",
            input_name,
            "\n".join(
                [
                    f"       {key}: {val}"
                    for key, val in component.get("filters", {}).items()
                ]
            ),
            "\n".join(
                [f"       {wildcard}" for wildcard in component.get("wildcards", [])]
            ),
        )
        return None
    try:
        path = itx.one(paths)
    except ValueError as err:
        msg = (
            f"More than one snakemake filename for {input_name}, taking the "
            f"first. To correct this, use the --filter_{input_name} option to "
            f"narrow the search. Found filenames: {paths}"
        )
        raise ConfigError(msg) from err

    if filters.has_empty_postfilter:
        return BidsComponent(
            name=input_name, path=path, zip_lists={key: [] for key in zip_lists}
        )

    return BidsComponent(name=input_name, path=path, zip_lists=zip_lists).filter(
        regex_search=True, **filters.post_exclusions
    )


def _parse_custom_path(
    input_path: Path | str,
    filters: UnifiedFilter,
) -> ZipList:
    """Glob wildcards from a custom path and apply filters.

    This replicates pybids path globbing for any custom path. Input path should have
    wildcards in braces as in "path/of/{wildcard_1}/{wildcard_2}_{wildcard_3}" Output
    will be arranged into a zip list of matches, list of matches, and Snakemake wildcard
    for each wildcard.

    Note that, currently, this will get confused if wildcard content matches
    non-wildcard content. For example, considering the path template above, the example
    "path/of/var1/variable_2_var3" would bug out because of the extra underscore.

    Parameters
    ----------
    input_path : str
        Path to be globbed
    regex_search : bool
        If True, use regex matching for filtering rather than simple equality
    **filters : str or list of str
        Values to keep. Each argument is the name of the entity to search

    Returns
    -------
    input_zip_list, input_list, input_wildcards
    """
    if not (wildcards := glob_wildcards(input_path)):
        _logger.warning("No wildcards defined in %s", input_path)

    # Log an error if no matches found
    if len(itx.first(wildcards.values())) == 0:
        _logger.error("No matching files for %s", input_path)
        return wildcards

    # Return the output values, running filtering on the zip_lists
    result = filter_list(
        wildcards,
        filters.without_bools,
        regex_search=False,
    )
    if not filters.post_exclusions:
        return result

    return filter_list(result, filters.post_exclusions, regex_search=True)


def _parse_bids_path(path: str, entities: Iterable[str]) -> tuple[str, dict[str, str]]:
    """Replace parameters in an bids path with the given wildcard {tags}.

    Parameters
    ----------
    path : str
        BIDS path
        (e.g. "root/sub-01/ses-01/sub-01_ses-01_T1w.nii.gz")
    wildcards : iterable of str
        BIDS entities to replace with wildcards. (e.g. "subject", "session", "suffix")

    Returns
    -------
    path : str
        Original path with the original entities replaced with wildcards.
        (e.g. "root/sub-{subject}/ses-{session}/sub-{subject}_ses-{session}_{suffix}")
    matches : iterable of (wildcard, value)
        The values matched with each wildcard
    """
    # If path is relative, we need to get a slash in front of it to ensure parsing works
    # correctly. So prepend "./" or ".\" and run function again, then strip before
    # returning
    if not os.path.isabs(path) and get_first_dir(path) != ".":
        path_, wildcard_values = _parse_bids_path(os.path.join(".", path), entities)
        return str(Path(path_)), wildcard_values

    entities = list(entities)

    matches = sorted(
        (
            (entity, match)
            for entity in map(BidsEntity, entities)
            for match in re.finditer(entity.regex, path)
        ),
        key=lambda match: match[1].start(2),
    )

    wildcard_values: dict[str, str] = {
        entity.wildcard: match.group(2) for entity, match in matches
    }
    if len(wildcard_values) != len(entities):
        unmatched = (
            set(map(BidsEntity, entities))
            .difference(match[0] for match in matches)
            .pop()
        )
        raise BidsParseError(path=path, entity=unmatched)

    num_matches = len(matches)
    new_path: list[str] = []
    for i in range(num_matches + 1):
        start = matches[i - 1][1].end(2) if i else 0
        end = len(path) if i == num_matches else matches[i][1].start(2)
        # Pybids technically allows `{` in the extension, so we escape it
        new_path.append(path[start:end].replace("{", "{{").replace("}", "}}"))
        if i < num_matches:
            new_path.append(f"{{{matches[i][0].wildcard}}}")
    return "".join(new_path), wildcard_values


def get_wildcard_constraints(image_types: InputsConfig) -> dict[str, str]:
    """Return a wildcard_constraints dict for use in snakemake.

    Contains constraints for all wildcards in the dynamically grabbed inputs.

    Parameters
    ----------
    image_types
        Component configuration dict

    Returns
    -------
        Dict containing wildcard constraints for all wildcards in the
        inputs, with typical bids naming constraints, ie letters and numbers
        ``[a-zA-Z0-9]+``.
    """
    bids_constraints = "[a-zA-Z0-9]+"
    return {
        entity: bids_constraints
        for imgtype in image_types
        for entity in image_types[imgtype]
    }
