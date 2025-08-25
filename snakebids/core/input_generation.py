# pyright: basic, reportPrivateUsage=false, reportUnknownMemberType=false, reportUnknownArgumentType=false, reportUnknownVariableType=false
"""Utilities for converting Snakemake apps to BIDS apps."""

from __future__ import annotations

import contextlib
import json
import logging
import os
import re
import warnings
from collections import defaultdict
from pathlib import Path
from typing import (
    Any,
    Callable,
    Iterable,
    Iterator,
    Literal,
    Mapping,
    overload,
)

import more_itertools as itx
from bids import BIDSLayout, BIDSLayoutIndexer

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
from snakebids.snakemake_compat import Snakemake
from snakebids.types import InputConfig, InputsConfig, ZipList
from snakebids.utils.snakemake_io import glob_wildcards
from snakebids.utils.utils import (
    DEPRECATION_FLAG,
    BidsEntity,
    BidsParseError,
    get_first_dir,
)

normalize_inputs_with_snakenull: Callable[..., Mapping[str, Any]] | None = None
try:
    from snakebids.snakenull import (
        normalize_inputs_with_snakenull as _normalize_inputs_with_snakenull,
    )

    normalize_inputs_with_snakenull = _normalize_inputs_with_snakenull
except (ImportError, ModuleNotFoundError):  # pragma: no cover
    # Leave normalize_inputs_with_snakenull as None
    pass

_logger = logging.getLogger(__name__)
_TOKEN_EQ_PLACEHOLDER = re.compile(r"^([a-zA-Z0-9]+)-\{\1\}$")


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
    snakenull: Mapping[str, Any] | None = ...,  # NEW in overload
) -> BidsDataset: ...


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
    snakenull: Mapping[str, Any] | None = ...,  # NEW in overload
) -> BidsDatasetDict: ...


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
    snakenull: Mapping[str, Any] | None = None,  # NEW in implementation
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
        snakenull_global=snakenull,
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
            stacklevel=3,
        )
    if pybids_reset_database is not None:
        warnings.warn(
            "The parameter `pybids_reset_database` in generate_inputs() is deprecated "
            "and will be removed in the next release. To reset the pybids database, "
            "use the `pybidsdb_reset` parameter instead.",
            stacklevel=3,
        )

    pybidsdb_dir = pybidsdb_dir or pybids_database_dir
    pybidsdb_reset = (
        pybidsdb_reset
        if pybidsdb_reset is not None
        else (pybids_reset_database if pybids_reset_database is not None else False)
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
            stacklevel=3,
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
                stacklevel=3,
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


def _is_local_relative(path: Path | str):
    """Test if a path location is local path relative to the current working directory.

    Parameter
    ---------
        path
            A UPath, Path, or str object to be checked

    Returns
    -------
    is_url : bool
        True if the path is relative
    """
    path_str = str(path)
    is_doubleslash_schemed = "://" in path_str
    return not is_doubleslash_schemed and not os.path.isabs(path_str)


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
        A boolean that determines whether to parse and index metadata

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


def _split_template(path_template: str) -> tuple[str, str]:
    """Return (dir, filename) split for a template path."""
    if "/" in path_template:
        head, tail = path_template.rsplit("/", 1)
        return head, tail
    return "", path_template


def _token_tags_in_filename(filename: str) -> set[str]:
    """Extract tags that appear as tokens 'tag-{tag}' in a filename."""
    tags: set[str] = set()
    for tok in filename.split("_"):
        m = _TOKEN_EQ_PLACEHOLDER.match(tok)
        if m:
            tags.add(m.group(1))
    return tags


def _drop_optional_tokens(filename: str, optional_tags: Iterable[str]) -> str:
    """Remove tokens 'tag-{tag}' for any tag in optional_tags from a filename."""
    opt = set(optional_tags)
    kept: list[str] = []
    for tok in filename.split("_"):
        m = _TOKEN_EQ_PLACEHOLDER.match(tok)
        if m and m.group(1) in opt:
            continue
        kept.append(tok)
    return "_".join(kept)


def _collapse_bids_templates(paths: set[str]) -> str | None:
    """Collapse multiple filename shapes by removing optional 'tag-{tag}' tokens.

    Optional = tokens that appear in some but not all templates. Only filename
    tokens are considered; directory segments like 'ses-{session}' are untouched.
    Returns a single unified template if all filenames become identical after
    removal; otherwise None.
    """
    if not paths or len(paths) == 1:
        return next(iter(paths), None)

    dir_files = [_split_template(p) for p in paths]
    dir_set = {d for d, _ in dir_files}
    if len(dir_set) != 1:
        # Different directory structures; don't attempt to collapse
        return None
    shared_dir = next(iter(dir_set))
    filenames = [f for _, f in dir_files]

    # Count presence of each 'tag-{tag}' token across filenames
    counts: dict[str, int] = {}
    for fn in filenames:
        for t in _token_tags_in_filename(fn):
            counts[t] = counts.get(t, 0) + 1

    n = len(filenames)
    optional_tags = {t for t, c in counts.items() if 0 < c < n}
    if not optional_tags:
        return None

    collapsed_files = [_drop_optional_tokens(fn, optional_tags) for fn in filenames]
    if len(set(collapsed_files)) == 1:
        collapsed_file = collapsed_files[0]
        return f"{shared_dir}/{collapsed_file}" if shared_dir else collapsed_file
    return None


def _snakenull_enabled_for_component(
    component_cfg: InputConfig, snakenull_global: Mapping[str, Any] | None
) -> bool:
    """Return True if snakenull is enabled for this component."""
    enabled = False
    if isinstance(snakenull_global, Mapping):
        val = snakenull_global.get("enabled")
        if isinstance(val, bool):
            enabled = val
    local = component_cfg.get("snakenull")
    if isinstance(local, Mapping):
        val = local.get("enabled")
        if isinstance(val, bool):
            enabled = val
    return enabled


def _maybe_normalize_component(
    name: str,
    comp: BidsComponent,
    inputs_config: InputsConfig,
    snakenull_global: Mapping[str, Any] | None,
) -> None:
    """Run snakenull normalizer for this component if enabled."""
    if normalize_inputs_with_snakenull is None:
        return
    if not _snakenull_enabled_for_component(inputs_config[name], snakenull_global):
        return
    cfg_for_this: dict[str, Any] = {"pybids_inputs": {name: inputs_config[name]}}
    if snakenull_global is not None:
        cfg_for_this["snakenull"] = dict(snakenull_global)
    normalize_inputs_with_snakenull({name: comp}, config=cfg_for_this)  # type: ignore[misc]


def _get_components(
    *,
    bids_layout: BIDSLayout | None,
    inputs_config: InputsConfig,
    postfilters: PostFilter,
    limit_to: Iterable[str] | None = None,
    snakenull_global: Mapping[str, Any] | None = None,
) -> Iterator[BidsComponent]:
    names = list(limit_to) if limit_to is not None else list(inputs_config.keys())

    for name in names:
        comp = _get_component(
            bids_layout=bids_layout,
            component=inputs_config[name],
            input_name=name,
            postfilters=postfilters,
            allow_template_collapse=_snakenull_enabled_for_component(
                inputs_config[name], snakenull_global
            ),
        )
        if comp is None:
            continue

        _maybe_normalize_component(name, comp, inputs_config, snakenull_global)
        yield comp


def _placeholders_in_template(path: str) -> list[str]:
    """Return placeholders that appear in the final path template (stable order)."""
    return list(dict.fromkeys(m.group(1) for m in re.finditer(r"{([^}]+)}", path)))


def _safe_parse_per_file(
    img, requested_wildcards: list[str], bids_layout
) -> tuple[str, dict[str, str]]:
    """Parse one file's path and return (template_path, per-file wildcard dict)."""
    parse_wildcards = [w for w in requested_wildcards if w in img.entities]
    _logger.debug("Wildcards %s found entities for %s", parse_wildcards, img.path)
    try:
        path, parsed_wildcards = _parse_bids_path(img.path, parse_wildcards)
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

    # parsed_wildcards uses wildcard names as keys, map to entity names
    per_file = {}
    for wc in requested_wildcards:
        wildcard_name = BidsEntity(wc).wildcard
        per_file[wc] = parsed_wildcards.get(wildcard_name, "")
    return path, per_file


def _build_zip_lists_from_parsed(
    parsed_per_file: list[dict[str, str]], placeholders: list[str]
) -> dict[str, list[str]]:
    """Rectangular zip_lists aligned with final template placeholders."""
    z: defaultdict[str, list[str]] = defaultdict(list)
    for rec in parsed_per_file:
        for wc in placeholders:
            # Find the entity name that corresponds to this wildcard name
            entity_name = None
            for entity_key in rec:
                if BidsEntity(entity_key).wildcard == wc:
                    entity_name = entity_key
                    break
            value = rec.get(entity_name, "") if entity_name else ""
            z[wc].append(value)
    return dict(z)


def _annotate_snakenull_component(
    comp, requested_wildcards: list[str], parsed_per_file: list[dict[str, str]]
) -> None:
    """Attach lightweight context for the inlined snakenull normalizer."""
    with contextlib.suppress(Exception):
        comp.snakenull_total_files = len(parsed_per_file)
    with contextlib.suppress(Exception):
        comp.requested_wildcards = requested_wildcards
    with contextlib.suppress(Exception):
        comp.snakenull_entities_per_file = parsed_per_file


def _select_template_path(
    paths: set[str], input_name: str, allow_collapse: bool
) -> str:
    """Return a single template path, optionally collapsing optional tokens."""
    if len(paths) == 1:
        return next(iter(paths))
    if allow_collapse:
        collapsed = _collapse_bids_templates(paths)
        if collapsed is not None:
            return collapsed
    msg = (
        f"Multiple path templates for one component. Use --filter_{input_name} "
        f"to narrow your search or --wildcards_{input_name} to make the template "
        "more generic.\n"
        f"\tcomponent = {input_name!r}\n"
        f"\tpath_templates = [\n\t\t" + ",\n\t\t".join(map(repr, paths)) + "\n\t]\n"
    ).expandtabs(4)
    raise ConfigError(msg)


def _make_component(
    input_name: str,
    path: str,
    zip_lists: dict[str, list[str]],
    filters: UnifiedFilter,
) -> BidsComponent:
    """Construct BidsComponent, respecting empty postfilters."""
    if filters.has_empty_postfilter:
        return BidsComponent(
            name=input_name, path=path, zip_lists={key: [] for key in zip_lists}
        )
    return BidsComponent(name=input_name, path=path, zip_lists=zip_lists).filter(
        regex_search=True, **filters.post_exclusions
    )


def _get_component(
    bids_layout: BIDSLayout | None,
    component: InputConfig,
    *,
    input_name: str,
    postfilters: PostFilter,
    allow_template_collapse: bool = False,
) -> BidsComponent | None:
    """Create component based on provided config."""
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

    try:
        matching_files = get_matching_files(bids_layout, filters)
    except FilterSpecError as err:
        raise err.get_config_error(input_name) from err

    requested_wildcards: list[str] = list(dict.fromkeys(component.get("wildcards", [])))

    # Parse each file; collect per-file dicts and candidate templates
    parsed_per_file: list[dict[str, str]] = []
    paths: set[str] = set()
    for img in matching_files:
        path, per_file = _safe_parse_per_file(img, requested_wildcards, bids_layout)
        parsed_per_file.append(per_file)
        paths.add(path)

    if not paths:
        _logger.warning(
            "No input files found for snakebids component %s:\n"
            "    filters:\n%s\n"
            "    wildcards:\n%s",
            input_name,
            "\n".join(
                f"       {key}: {val}"
                for key, val in component.get("filters", {}).items()
            ),
            "\n".join(f"       {wc}" for wc in requested_wildcards),
        )
        return None

    # Decide on a single template (collapse only if allowed)
    path = _select_template_path(paths, input_name, allow_template_collapse)

    # Build rectangular zip_lists only for placeholders in the final path
    placeholders = _placeholders_in_template(path)
    zip_lists = _build_zip_lists_from_parsed(parsed_per_file, placeholders)

    # Construct component and annotate context for snakenull
    comp = _make_component(input_name, path, zip_lists, filters)
    _annotate_snakenull_component(comp, requested_wildcards, parsed_per_file)
    return comp


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
    if _is_local_relative(path) and get_first_dir(path) != ".":
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
