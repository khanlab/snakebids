"""Utilities for converting Snakemake apps to BIDS apps."""
from __future__ import annotations

import json
import logging
import re
from collections import defaultdict
from pathlib import Path
from typing import Any, Generator, Iterable, Optional, Sequence, overload

import more_itertools as itx
from bids import BIDSLayout, BIDSLayoutIndexer
from bids.layout import BIDSFile, Query
from snakemake.script import Snakemake
from typing_extensions import Literal

from snakebids.core.datasets import BidsComponent, BidsDataset, BidsDatasetDict
from snakebids.core.filtering import filter_list
from snakebids.exceptions import ConfigError, PybidsError
from snakebids.types import InputsConfig, ZipList
from snakebids.utils.snakemake_io import glob_wildcards
from snakebids.utils.utils import BidsEntity, BidsParseError, MultiSelectDict, surround

_logger = logging.getLogger(__name__)


@overload
def generate_inputs(  # noqa: PLR0913
    bids_dir: Path | str,
    pybids_inputs: InputsConfig,
    pybids_database_dir: Path | str | None = ...,
    pybids_reset_database: bool = ...,
    derivatives: bool | Path | str = ...,
    pybids_config: str | None = ...,
    limit_to: Iterable[str] | None = ...,
    participant_label: Iterable[str] | str | None = ...,
    exclude_participant_label: Iterable[str] | str | None = ...,
    use_bids_inputs: Literal[True] | None = ...,
) -> BidsDataset:
    ...


@overload
def generate_inputs(  # noqa: PLR0913
    bids_dir: Path | str,
    pybids_inputs: InputsConfig,
    pybids_database_dir: Path | str | None = ...,
    pybids_reset_database: bool = ...,
    derivatives: bool | Path | str = ...,
    pybids_config: str | None = ...,
    limit_to: Iterable[str] | None = ...,
    participant_label: Iterable[str] | str | None = ...,
    exclude_participant_label: Iterable[str] | str | None = ...,
    use_bids_inputs: Literal[False] = ...,
) -> BidsDatasetDict:
    ...


def generate_inputs(  # noqa: PLR0913
    bids_dir: Path | str,
    pybids_inputs: InputsConfig,
    pybids_database_dir: Path | str | None = None,
    pybids_reset_database: bool = False,
    derivatives: bool | Path | str = False,
    pybids_config: str | None = None,
    limit_to: Iterable[str] | None = None,
    participant_label: Iterable[str] | str | None = None,
    exclude_participant_label: Iterable[str] | str | None = None,
    use_bids_inputs: bool | None = None,
) -> BidsDataset | BidsDatasetDict:
    """Dynamically generate snakemake inputs using pybids_inputs

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

    pybids_database_dir
        Path to database directory. If None is provided, database
        is not used

    pybids_reset_database
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

    subject_filter, regex_search = _generate_filters(
        participant_label, exclude_participant_label
    )

    # Generates a BIDSLayout
    layout = (
        _gen_bids_layout(
            bids_dir=bids_dir,
            derivatives=derivatives,
            pybids_config=pybids_config,
            pybids_database_dir=pybids_database_dir,
            pybids_reset_database=pybids_reset_database,
        )
        if not _all_custom_paths(pybids_inputs)
        else None
    )

    filters = {"subject": subject_filter} if subject_filter else {}
    bids_inputs = _get_lists_from_bids(
        bids_layout=layout,
        pybids_inputs=pybids_inputs,
        limit_to=limit_to,
        regex_search=regex_search,
        **(filters),
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
    except ValueError as err:
        raise ConfigError(
            "Multiple components found with the same name: "
            + ", ".join(surround(err.args[0], "'"))
        ) from err

    if use_bids_inputs:
        return dataset
    return dataset.as_dict


def _all_custom_paths(config: InputsConfig):
    """Check that all input components have a custom path"""
    return all(comp.get("custom_path") for comp in config.values())


def _gen_bids_layout(
    bids_dir: Path | str,
    derivatives: Path | str | bool,
    pybids_database_dir: Path | str | None,
    pybids_reset_database: bool,
    pybids_config: Path | str | None = None,
) -> BIDSLayout:
    """Create (or reindex) the BIDSLayout if one doesn't exist,
    which is only saved if a database directory path is provided

     Parameters
    ----------
    bids_dir
        Path to bids directory

    derivatives
        A boolean (or path(s) to derivatives datasets) that
        determines whether snakebids will search in the
        derivatives subdirectory of the input dataset.

    pybids_database_dir
        Path to database directory. If None is provided, database
        is not used

    pybids_reset_database
        A boolean that determines whether to reset / overwrite
        existing database.

    Returns
    -------
    layout : BIDSLayout
        Layout from pybids for accessing the BIDS dataset to grab paths
    """

    # Check for database_dir
    # If blank, assume db not to be used
    if not pybids_database_dir:
        pybids_database_dir = None
    # Otherwise check for relative path and update
    elif not Path(pybids_database_dir).is_absolute():
        pybids_database_dir = None
        _logger.warning("Absolute path must be provided, database will not be used")

    return BIDSLayout(
        str(bids_dir),
        derivatives=derivatives,
        validate=False,
        config=pybids_config,
        database_path=pybids_database_dir,
        reset_database=pybids_reset_database,
        indexer=BIDSLayoutIndexer(validate=False, index_metadata=False),
    )


def write_derivative_json(snakemake: Snakemake, **kwargs: dict[str, Any]) -> None:
    """Snakemake function to read a json file, and write to a new one,
    adding BIDS derivatives fields for Sources and Parameters.

    Parameters
    ----------
    snakemake : struct Snakemake
        structure passed to snakemake python scripts, containing input,
        output, params, config ...
        This function requires input.json and output.json to be defined, as
        it will read and write json files
    """

    with open(snakemake.input.json, "r", encoding="utf-8") as input_json:
        sidecar = json.load(input_json)

    sidecar.update(
        {
            "Sources": [snakemake.input],
            "Parameters": snakemake.params,
            **kwargs,
        }
    )

    with open(snakemake.output.json, "w", encoding="utf-8") as outfile:
        json.dump(sidecar, outfile, indent=4)


def _generate_filters(
    include: Iterable[str] | str | None = None,
    exclude: Iterable[str] | str | None = None,
) -> tuple[list[str], bool]:
    """Generate Pybids filter based on inclusion or exclusion criteria

    Converts either a list of values to include or exclude in a list of Pybids
    compatible filters. Unlike inclusion values, exclusion requires regex filtering. The
    necessity for regex will be indicated by the boolean value of the second returned
    item: True if regex is needed, False otherwise. Raises an exception if both include
    and exclude are stipulated

    Parameters
    ----------
    include : list of str or str, optional
        Values to include, values not found in this list will be excluded, by default
        None
    exclude : list of str or str, optional
        Values to exclude, only values not found in this list will be included, by
        default None

    Returns
    -------
    list of str, bool
        Two values: the first, a list of pybids compatible filters; the second, a
        boolean indicating whether regex_search must be enabled in pybids

    Raises
    ------
    ValueError Raised of both include and exclude values are stipulated.
    """
    if include is not None and exclude is not None:
        raise ValueError(
            "Cannot define both participant_label and "
            "exclude_participant_label at the same time"
        )

    # add participant_label or exclude_participant_label to search terms (if
    # defined)
    # we make the item key in search_terms a list so we can have both
    # include and exclude defined
    if include is not None:
        return [*itx.always_iterable(include)], False

    if exclude is not None:
        # if multiple items to exclude, combine with with item1|item2|...
        exclude_string = "|".join(
            re.escape(label) for label in itx.always_iterable(exclude)
        )
        # regex to exclude subjects
        return [f"^((?!({exclude_string})$).*)$"], True
    return [], False


def _parse_custom_path(
    input_path: Path | str,
    regex_search: bool = False,
    **filters: list[str] | str,
) -> ZipList:
    """Glob wildcards from a custom path and apply filters

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
    wildcards = glob_wildcards(input_path)

    if not wildcards:
        _logger.warning("No wildcards defined in %s", input_path)

    # Log an error if no matches found
    if len(next(iter(wildcards.values()))) == 0:
        _logger.error("No matching files for %s", input_path)
        return wildcards

    # Return the output values, running filtering on the zip_lists
    return filter_list(wildcards, filters, regex_search=regex_search)


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

    wildcard_values: dict[str, str] = {}

    for entity in map(BidsEntity, entities):
        # Iterate over wildcards, slowly updating the path as each entity is replaced

        wildcard = entity.wildcard
        match = re.search(entity.regex, path)
        if not match or not match.group(2):
            raise BidsParseError(path=path, entity=entity)

        # overwrite path one wildcard at a time
        path = re.sub(entity.regex, rf"\1{{{wildcard}}}\3", path)
        value = match.group(2)
        wildcard_values[wildcard] = value

    return path, wildcard_values


def _get_lists_from_bids(
    bids_layout: Optional[BIDSLayout],
    pybids_inputs: InputsConfig,
    *,
    limit_to: Iterable[str] | None = None,
    regex_search: bool = False,
    **filters: str | Sequence[str],
) -> Generator[BidsComponent, None, None]:
    """Grabs files using pybids and creates snakemake-friendly lists

    Parameters
    ----------
    bids_layout : BIDSLayout
        Layout from pybids for accessing the BIDS dataset to grab paths

    pybids_inputs : dict
        Dictionary indexed by modality name, specifying the filters and
        wildcards for each pybids input.

    limit_to : list, optional
        List of inputs to skip, this used by snakebids to exclude modalities
        based on cmd-line args

    filters : dict of str -> str or list of str, optional
        Pybids filters to apply globally to all inputs.

    Yields
    ------
    BidsComponent:
        One BidsComponent is yielded for each modality described by ``pybids_inputs``.
    """
    for input_name in limit_to or list(pybids_inputs):
        _logger.debug("Grabbing inputs for %s...", input_name)
        component = pybids_inputs[input_name]

        if "custom_path" in component:
            # a custom path was specified for this input, skip pybids:
            # get input_wildcards by parsing path for {} entries (using a set
            # to get unique only)
            # get zip_lists by using glob_wildcards (but need to modify
            # to deal with multiple wildcards

            path = component["custom_path"]
            zip_lists = _parse_custom_path(
                path,
                regex_search=regex_search,
                **pybids_inputs[input_name].get("filters", {}),
                **filters,
            )
            yield BidsComponent(input_name, path, MultiSelectDict(zip_lists))
            continue

        if bids_layout is None:
            _logger.warning(
                "No valid bids dir given, but %s does not have a custom_path specified "
                "and will be skipped.",
                input_name,
            )
            continue

        zip_lists: dict[str, list[str]] = defaultdict(list)
        paths: set[str] = set()
        pybids_filters = {
            key: Query.ANY if val is True else Query.NONE if val is False else val
            for key, val in component.get("filters", {}).items()
        }
        try:
            matching_files: Iterable[BIDSFile] = bids_layout.get(
                regex_search=regex_search, **pybids_filters, **filters
            )
        except AttributeError as err:
            raise PybidsError(
                "Pybids has encountered a problem that Snakebids cannot handle. This "
                "may indicate a missing or invalid dataset_description.json for this "
                "dataset."
            ) from err
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
                raise ConfigError(
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
                ) from err

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
                    [
                        f"       {wildcard}"
                        for wildcard in component.get("wildcards", [])
                    ]
                ),
            )
            continue
        try:
            path = itx.one(paths)
        except ValueError:
            raise ConfigError(
                f"More than one snakemake filename for {input_name}, taking the "
                f"first. To correct this, use the --filter_{input_name} option to "
                f"narrow the search. Found filenames: {paths}"
            )

        yield BidsComponent(input_name, path, MultiSelectDict(zip_lists))


def get_wildcard_constraints(image_types: InputsConfig) -> dict[str, str]:
    """Return a wildcard_constraints dict for snakemake to use, containing
    all the wildcards that are in the dynamically grabbed inputs

    Parameters
    ----------
    image_types : dict

    Returns
    -------
        Dict containing wildcard constraints for all wildcards in the
        inputs, with typical bids naming constraints, ie letters and numbers
        ``[a-zA-Z0-9]+``.
    """
    bids_constraints = "[a-zA-Z0-9]+"
    return {
        entity: bids_constraints
        for imgtype in image_types.keys()
        for entity in image_types[imgtype].keys()
    }
