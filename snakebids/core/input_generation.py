import json
import logging
import operator as op
import re
from pathlib import Path
from typing import Dict, List, Tuple, Union

import more_itertools as itx
from bids import BIDSLayout, BIDSLayoutIndexer

from snakebids.core.filtering import filter_list
from snakebids.utils.snakemake_io import glob_wildcards
from snakebids.utils.utils import get_match_search_func, read_bids_tags

_logger = logging.getLogger(__name__)


# pylint: disable=too-many-arguments,too-many-locals
def generate_inputs(
    bids_dir,
    pybids_inputs,
    pybids_database_dir=None,
    pybids_reset_database=False,
    derivatives=False,
    pybids_config=None,
    limit_to=None,
    participant_label=None,
    exclude_participant_label=None,
):
    """Dynamically generate snakemake inputs using pybids_inputs dict, and
    pybids to parse the bids dataset.

    Parameters
    ----------
    bids_dir : str
        Path to bids directory

    pybids_inputs : dict
        Configuration for bids inputs, with keys as the names (``str``)

        Nested `dicts` with the following required keys:

        * ``"filters"``: Dictionary containing keyword arguments that will
          be passed to pybids ``get()``.

        * ``"wildcards"``: List of (str) bids tags to include as wildcards in
          snakemake. At minimum this should usually include
          ``['subject','session']``, plus any other wildcards that you may
          want to make use of in your snakemake workflow, or want to retain
          in the output paths. Any wildcards in this list that are not in the
          filename will just be ignored.

    Returns
    -------
    dict:
        The dict returned by this functions contains seven items. Each of
        the following four items is a dict containing one item for each
        modality described by ``pybids_inputs``.

        * ``"input_path"``: String with a wildcard-filled path that matches
          the images for this modality.

        * ``"input_zip_lists"``: Dictionary where each key is a wildcard
          entity and each value is a list of the values found for that
          entity. Each of these lists has length equal to the number of
          images matched for this modality, so they can be zipped together to
          get a list of the wildcard values for each image.

        * ``"input_lists"``: Dictionary where each key is a wildcard entity
          and each value is a list of the unique values found for that
          entity. These lists may not be the same length.

        * ``"input_wildcards"``: Dictionary where each key is the name of a
          wildcard entity, and each value is the Snakemake wildcard used for
          that entity.

        Then there are three more top-level entries in the dictionary:

        * ``"subjects"``: A list of the subjects in the dataset.

        * ``"sessions"``: A list of the sessions in the dataset.

        * ``"subj_wildcards"``: The subject and session wildcards applicable
          to this dataset. ``{"subject": "{subject}"}`` if there is only one
          session, ``{"subject": "{subject}", "session": "{session}"}`` if
          there are multiple sessions.

    Notes
    -----
    As an example, consider the following BIDS dataset::

        bids-example/
        ├── dataset_description.json
        ├── participants.tsv
        ├── README
        └── sub-control01
            ├── anat
            │   ├── sub-control01_T1w.json
            │   ├── sub-control01_T1w.nii.gz
            │   ├── sub-control01_T2w.json
            │   └── sub-control01_T2w.nii.gz
            ├── dwi
            │   ├── sub-control01_dwi.bval
            │   ├── sub-control01_dwi.bvec
            │   └── sub-control01_dwi.nii.gz
            ├── fmap
            │   ├── sub-control01_magnitude1.nii.gz
            │   ├── sub-control01_phasediff.json
            │   ├── sub-control01_phasediff.nii.gz
            │   └── sub-control01_scans.tsv
            └── func
                ├── sub-control01_task-nback_bold.json
                ├── sub-control01_task-nback_bold.nii.gz
                ├── sub-control01_task-nback_events.tsv
                ├── sub-control01_task-nback_physio.json
                ├── sub-control01_task-nback_physio.tsv.gz
                ├── sub-control01_task-nback_sbref.nii.gz
                ├── sub-control01_task-rest_bold.json
                ├── sub-control01_task-rest_bold.nii.gz
                ├── sub-control01_task-rest_physio.json
                └── sub-control01_task-rest_physio.tsv.gz

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

    Then ``generate_inputs(bids_dir, pybids_input)`` would return the
    following dictionary::

        {
            "input_path": {
                "bold": "bids-example/sub-{subject}/func/sub-{subject}_task-{task}_bold.nii.gz"
            },
            "input_zip_lists": {
                "bold": {
                    "subject": ["control01", "control01"],
                    "task": ["nback", "rest"]
                }
            },
            "input_lists": {
                "bold": {
                    "subject": ["control01"],
                    "task": ["nback", "rest"]
                }
            },
            "input_wildcards": {
                "bold": {
                    "subject": "{subject}",
                    "task": "{task}"
                }
            },
            "subjects": ["subject01"],
            "sessions": [],
            "subj_wildcards": {"subject": "{subject}"}
        }
    """  # noqa

    subject_filter, regex_search = _generate_filters(
        participant_label, exclude_participant_label
    )

    # Generates a BIDSLayout
    layout = _gen_bids_layout(
        bids_dir=bids_dir,
        derivatives=derivatives,
        pybids_config=pybids_config,
        pybids_database_dir=pybids_database_dir,
        pybids_reset_database=pybids_reset_database,
    )

    # this will populate input_path, input_lists, input_zip_lists, and
    # input_wildcards

    filters = {"subject": subject_filter} if subject_filter else {}
    inputs_config_dict = _get_lists_from_bids(
        bids_layout=layout,
        pybids_inputs=pybids_inputs,
        limit_to=limit_to,
        regex_search=regex_search,
        **filters,
    )

    if layout is None:
        # if no layout, then use subjects/sessions from --path vars
        subjects = []
        sessions = []
        for input_type in inputs_config_dict["input_lists"]:

            subj_set = set(inputs_config_dict["input_lists"][input_type]["subject"])

            # filter the list of subjects with participant_label
            if participant_label is not None:
                subj_set = set.intersection(subj_set, set(participant_label))

            # TODO: need to also remove subjects based on exclude_participant_label

            # replace with filtered list
            inputs_config_dict["input_lists"][input_type]["subject"] = list(subj_set)

            # add to set of subjects from all input_types
            subjects.append(subj_set)

            if "session" in (inputs_config_dict["input_lists"][input_type].keys()):
                sessions.append(
                    {inputs_config_dict["input_lists"][input_type]["session"]}
                )
            else:
                sessions.append(set([]))

        # take set intersection of all input types
        inputs_config_dict["subjects"] = list(set.intersection(*subjects))
        inputs_config_dict["sessions"] = list(set.intersection(*sessions))

    else:
        # populate subjects, sessions and subj_wildcards in the config
        inputs_config_dict["subjects"] = layout.get_subjects(**filters)
        inputs_config_dict["sessions"] = layout.get_sessions(**filters)

    if len(inputs_config_dict["sessions"]) == 0:
        inputs_config_dict["subj_wildcards"] = {"subject": "{subject}"}
    else:
        inputs_config_dict["subj_wildcards"] = {
            "subject": "{subject}",
            "session": "{session}",
        }

    return inputs_config_dict


def _gen_bids_layout(
    bids_dir,
    derivatives,
    pybids_database_dir,
    pybids_reset_database,
    pybids_config=None,
):
    """Create (or reindex) the BIDSLayout if one doesn't exist,
    which is only saved if a database directory path is provided
    """

    # Set db dir to None (otherwise saves to parent dir)
    if Path(bids_dir):
        # Check for database_dir
        # If blank, assume db not to be used
        if pybids_database_dir == "":
            pybids_database_dir = None
        # Otherwise check for relative path and update
        elif (
            pybids_database_dir is not None
            and not Path(pybids_database_dir).is_absolute()
        ):
            pybids_database_dir = None
            _logger.warning(
                "Absolute path must be provided, database will " "not be used"
            )

        layout = BIDSLayout(
            bids_dir,
            derivatives=derivatives,
            validate=False,
            config=pybids_config,
            database_path=pybids_database_dir,
            reset_database=pybids_reset_database,
            indexer=BIDSLayoutIndexer(validate=False, index_metadata=False),
        )
    else:
        _logger.info(
            "bids_dir does not exist, skipping PyBIDS and using "
            "custom file paths only"
        )
        layout = None

    return layout


def write_derivative_json(snakemake, **kwargs):
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
    include: Union[List[str], str] = None, exclude: Union[List[str], str] = None
) -> Tuple[List[str], bool]:
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
        if isinstance(exclude, list):
            exclude_string = "|".join(re.escape(label) for label in exclude)
        # if not, then string is the label itself
        else:
            exclude_string = re.escape(exclude)
        # regex to exclude subjects
        return [f"^((?!({exclude_string})$).*)$"], True
    return [], False


def _parse_custom_path(
    input_path: Union[Path, str],
    regex_search: bool = False,
    **filters: Union[List[str], str],
):
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
    wildcard_names = list(wildcards._fields)

    if len(wildcard_names) == 0:
        _logger.warning("No wildcards defined in %s", input_path)

    # Initialize output values
    input_wildcards: Dict[str, str] = {}
    input_zip_lists: Dict[str, Tuple[str]] = {}
    input_lists: Dict[str, List[str]] = {}

    # Loop through every wildcard name
    for i, wildcard in enumerate(wildcard_names):
        # Check if this wildcard needs to be filtered
        if wildcard not in filters:
            # If not, the match_func will always return True
            match_func = lambda _: True  # noqa: E731
        elif regex_search:
            # Otherwise, we use regex matching for regex_search
            match_func = get_match_search_func(
                itx.always_iterable(filters[wildcard]), re.match
            )
        else:
            # And a simple equality operator for everything else
            match_func = get_match_search_func(
                itx.always_iterable(filters[wildcard]), op.eq
            )

        # Add the wildcard item to each output value, using filtering for input_lists
        input_zip_lists[wildcard] = wildcards[i]
        input_lists[wildcard] = [*filter(match_func, set(wildcards[i]))]
        input_wildcards[wildcard] = f"{{{wildcard}}}"

        # Log an error if no matches found
        # TODO: This will fail to detect filtering correctly as, up till now, it has
        #       only been performed on input_lists
        if len(wildcards[i]) == 0:
            _logger.error("No matching files for %s", input_path)

    # Return the output values, running filtering on the input_zip_lists
    return (
        filter_list(input_zip_lists, filters, regex_search=regex_search),
        input_lists,
        input_wildcards,
    )


def _process_layout_wildcard(path, wildcard_name):
    """Convert an absolute BIDS path to the same path with the given tag
    replaced by a wildcard.

    Parameters
    ----------
    path : str
        Absolute BIDS path
        (e.g. "root/sub-01/ses-01/sub-01_ses-01_T1w.nii.gz")
    wildcard_name : str
        BIDS entity to replace with a wildcard. (e.g. "subject")

    Returns
    -------
    path : str
        Original path with the original entity replaced with a wildcard.
        (e.g. "root/sub-{subject}/ses-01/sub-{subject}_ses-01_T1w.nii.gz")
    match : str
        Matched BIDS entity in the input path (e.g. "01")
    out_name : str
        Name of the applied wildcard (e.g. "subject")
    """
    bids_tags = read_bids_tags()
    tag = bids_tags[wildcard_name] if wildcard_name in bids_tags else wildcard_name

    # this changes e.g. sub-001 to sub-{subject} in the path
    # (so snakemake can use the wildcards)
    # HACK FIX FOR acq vs acquisition etc -- should
    # eventually update the bids() function to also use
    # bids_tags.json, where e.g. acquisition -> acq is
    # defined.. -- then, can use wildcard_name instead
    # of out_name..
    if wildcard_name not in ["subject", "session"]:
        out_name = tag
    else:
        out_name = wildcard_name

    if wildcard_name == "suffix":
        # capture suffix
        matching_pattern = ".*_([a-zA-Z0-9]+).*$"
        # capture before and after suffix
        replace_pattern = "(.*_)[a-zA-Z0-9]+(.*)$"
        # replace with before, {suffix}, after
        replace = f"\\1{{{out_name}}}\\2"
        match = re.search(matching_pattern, path)
        path = re.sub(replace_pattern, replace, path)

    else:
        pattern = f"{tag}-([a-zA-Z0-9]+)"
        replace = f"{tag}-{{{out_name}}}"

        match = re.search(pattern, path)
        path = re.sub(pattern, replace, path)

    # update the path with the {wildcards} -- uses the
    # value from the string (not from the pybids
    # entities), since that has issues with integer
    # formatting (e.g. for run=01)

    return path, match[1], out_name


def _get_lists_from_bids(bids_layout, pybids_inputs, limit_to=None, **filters):
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

    filters : dict, optional
        Pybids filters to apply globally to all inputs.

    Returns
    -------
    dict of dicts:
        There is one entry in the dictionary for each modality described by
        ``pybids_inputs``. Each of these entries contains the following items:

        * ``"input_path"``: String with a wildcard-filled path that matches
          the images for this modality.
        * ``"input_zip_lists"``: Dictionary where each key is a wildcard
          entity and each value is a list of the values found for that
          entity. Each of these lists has length equal to the number of
          images matched for this modality, so they can be zipped together to
          get a list of the wildcard values for each image.
        * ``"input_lists"``: Dictionary where each key is a wildcard entity
          and each value is a list of the unique values found for that
          entity. These lists may not be the same length.
        * ``"input_wildcards"``: Dictionary where each key is the name of a
          wildcard entity, and each value is the Snakemake wildcard used for
          that entity.
    """

    out_dict = dict(
        {
            "input_path": {},
            "input_zip_lists": {},
            "input_lists": {},
            "input_wildcards": {},
        }
    )

    if limit_to is None:
        limit_to = pybids_inputs.keys()

    for input_name in limit_to:
        _logger.debug("Grabbing inputs for %s...", input_name)
        input_path = ""
        input_wildcards = {}
        input_zip_lists = {}
        input_lists = {}
        if "custom_path" in pybids_inputs[input_name].keys():
            # a custom path was specified for this input, skip pybids:
            # get input_wildcards by parsing path for {} entries (using a set
            # to get unique only)
            # get input_zip_lists by using glob_wildcards (but need to modify
            # to deal with multiple wildcards

            input_path = pybids_inputs[input_name]["custom_path"]
            (
                input_zip_lists,
                input_lists,
                input_wildcards,
            ) = _parse_custom_path(input_path)
        else:
            paths = set()
            for img in bids_layout.get(
                **pybids_inputs[input_name].get("filters", {}), **filters
            ):
                input_path = img.path
                for wildcard_name in pybids_inputs[input_name].get("wildcards", []):

                    if wildcard_name not in img.get_entities():
                        continue
                    _logger.debug(
                        "Wildcard %s found entities for %s",
                        wildcard_name,
                        img.path,
                    )
                    (
                        input_path,
                        input_list,
                        out_name,
                    ) = _process_layout_wildcard(input_path, wildcard_name)
                    if out_name not in input_zip_lists:
                        input_zip_lists[out_name] = []
                        input_lists[out_name] = set()
                        input_wildcards[out_name] = {}
                    input_zip_lists[out_name].append(input_list)
                    input_lists[out_name].add(input_list)
                    input_wildcards[out_name] = f"{{{out_name}}}"
                paths.add(input_path)

            # now, check to see if unique
            if len(paths) == 0:
                _logger.warning("No images found for %s", input_name)
                continue
            if len(paths) > 1:
                _logger.warning(
                    "More than one snakemake filename for %s, taking the "
                    "first. To correct this, use the --filter_%s option to "
                    "narrow the search. Found filenames: %s",
                    input_name,
                    input_name,
                    paths,
                )

            input_path = list(paths)[0]

            # convert sets to lists
            input_lists = {key: list(val) for key, val in input_lists.items()}

        out_dict["input_path"][input_name] = input_path
        out_dict["input_zip_lists"][input_name] = input_zip_lists
        out_dict["input_lists"][input_name] = input_lists
        out_dict["input_wildcards"][input_name] = input_wildcards

    return out_dict


def get_wildcard_constraints(image_types):
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
