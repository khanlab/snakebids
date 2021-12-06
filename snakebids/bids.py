"""Utilities for converting Snakemake apps to BIDS apps."""

import json
import logging
import os
import re
from collections import OrderedDict
from os.path import join
from pathlib import Path
from typing import Optional

from bids import BIDSLayout, BIDSLayoutIndexer

from snakebids.utils.snakemake_io import glob_wildcards

_logger = logging.getLogger(__name__)


# pylint: disable=too-many-arguments
def bids(
    root=None,
    datatype=None,
    prefix=None,
    suffix=None,
    subject=None,
    session=None,
    include_subject_dir=True,
    include_session_dir=True,
    **entities,
):
    """Helper function for generating bids paths for snakemake workflows.

    File path is of the form::

        [root]/[sub-{subject}]/[ses-{session]/
        [prefix]_[sub-{subject}]_[ses-{session}]_[{key}-{val}_ ... ]_[suffix]

    Parameters
    ----------
    root : str or Path, default=None
        root folder to include in the path (e.g. 'results')
    datatype : str, default=None
        folder to include after sub-/ses- (e.g. anat, dwi )
    prefix : str, default=None
        string to prepend to the file name (typically not defined, unless you
        want tpl-{tpl}, or a datatype)
    suffix : str, default=None
        bids suffix including extension (e.g. 'T1w.nii.gz')
    subject : str, default=None
        subject to use, for folder and filename
    session : str, default=None
        session to use, for folder and filename
    include_subject_dir : bool, default=True
        whether to include the sub-{subject} folder if subject defined
        (default: True)
    include_session_dir : bool, default=True
        whether to include the ses-{session} folder if session defined
        (default: True)
    **entities : dict, optional
        dictionary of bids entities (e.g. space=T1w for space-T1w)

    Returns
    -------
    str
        bids-like file path

    Examples
    --------

    Below is a rule using bids naming for input and output::

        rule proc_img:
            input: 'sub-{subject}_T1w.nii.gz'
            output: 'sub-{subject}_space-snsx32_desc-preproc_T1w.nii.gz'

    With bids() you can instead use::

         rule proc_img:
            input: bids(subject='{subject}',suffix='T1w.nii.gz')
            output: bids(
                subject='{subject}',
                space='snsx32',
                desc='preproc',
                suffix='T1w.nii.gz'
            )

    Note that here we are not actually using "functions as inputs" in
    snakemake, which would require a function definition with wildcards as
    the argument, and restrict to input/params, but bids() is being used
    simply to return a string.

    Also note that space, desc and suffix are NOT wildcards here, only
    {subject} is. This makes it easy to combine wildcards and non-wildcards
    with bids-like naming.

    However, you can still use bids() in a lambda function. This is
    especially useful if your wildcards are named the same as bids entities
    (e.g. {subject}, {session}, {task} etc..)::

        rule proc_img:
            input: lambda wildcards: bids(**wildcards,suffix='T1w.nii.gz')
            output: bids(
                subject='{subject}',
                space='snsx32',
                desc='preproc',
                suffix='T1w.nii.gz'
            )

    Or another example where you may have many bids-like wildcards used in
    your workflow::

        rule denoise_func:
            input: lambda wildcards: bids(**wildcards, suffix='bold.nii.gz')
            output: bids(
                subject='{subject}',
                session='{session}',
                task='{task}',
                acq='{acq}',
                desc='denoise',
                suffix='bold.nii.gz'
            )

    In this example, all the wildcards will be determined from the output
    and passed on to bids() for inputs. The output filename will have a
    'desc-denoise' flag added to it.

    Also note that even if you supply entities in a different order, the
    entities will be ordered based on the OrderedDict defined here.
    If entities not known are provided, they will be just be placed
    at the end (before the suffix), in the order you provide them in.

    Notes
    -----

    * For maximum flexibility all arguments are optional (if none are
      specified, will return empty string)

    * Some code adapted from mne-bids, specifically
      https://mne.tools/mne-bids/stable/_modules/mne_bids/utils.html
    """

    # replace underscores in keys (needed to that users can use reserved
    # keywords by appending a _)
    entities = {k.replace("_", ""): v for k, v in entities.items()}

    # strict ordering of bids entities is specified here:
    order: OrderedDict[str, Optional[str]] = OrderedDict(
        [
            ("task", None),
            ("acq", None),
            ("ce", None),
            ("rec", None),
            ("dir", None),
            ("run", None),
            ("mod", None),
            ("echo", None),
            ("hemi", None),
            ("space", None),
            ("res", None),
            ("den", None),
            ("label", None),
            ("desc", None),
        ]
    )

    # Now add in entities (this preserves ordering above)
    for key, val in entities.items():
        order[key] = val

    # initialize lists for filename and folder
    # will append to these, then '_'.join() os.path.join() respectively
    filename = []
    folder = []

    # root directory
    if isinstance(root, str):
        folder.append(root)
    elif isinstance(root, Path):
        folder.append(str(root.resolve()))

    # if prefix is defined, put it before other anything else
    if isinstance(prefix, str):
        filename.append(prefix)

    # if subject defined then append to file and folder
    if isinstance(subject, str):
        if include_subject_dir is True:
            folder.append(f"sub-{subject}")
        filename.append(f"sub-{subject}")

    # if session defined then append to file and folder
    if isinstance(session, str):
        if include_session_dir is True:
            folder.append(f"ses-{session}")
        filename.append(f"ses-{session}")

    if isinstance(datatype, str):
        folder.append(datatype)

    # add the entities
    filename += [f"{key}-{val}" for key, val in order.items() if val is not None]

    # if suffix is defined, append it
    if isinstance(suffix, str):
        filename.append(suffix)

    if len(filename) == 0:
        return ""

    # now, join up the lists:
    filename = "_".join(filename)

    if len(folder) > 0:
        filename = join(*folder, filename)

    return filename


def print_boilerplate():
    """Function to print out boilerplate to add to Snakefile. (not used
    anywhere yet)"""

    print(
        """
# ---- begin snakebids boilerplate ------------------------------------------

import snakebids
from snakebids import bids

configfile: 'config/snakebids.yml'

#writes inputs_config.yml and updates config dict
config.update(snakebids.generate_inputs(bids_dir=config['bids_dir'],
                            pybids_inputs=config['pybids_inputs'],
                            derivatives=config['derivatives'],
                            participant_label=config['participant_label'],
                            exclude_participant_label=config['exclude_participant_label']))


#this adds constraints to the bids naming
wildcard_constraints:  **snakebids.get_wildcard_constraints(
    config['pybids_inputs']
)

# ---- end snakebids boilerplate --------------------------------------------
"""
    )


def filter_list(zip_list, wildcards, return_indices_only=False):
    """This function is used when you are expanding over some subset of the
    wildcards i.e. if your output file doesn't contain all the wildcards in
    input_wildcards

    Parameters
    ----------
    zip_list : dict
        generated zip lists dict from config file to filter

    wildcards : dict
        wildcard values to filter the zip lists

    return_indices_only : bool, default=False
        return the indices of the matching wildcards

    Returns
    -------
    dict
        zip list with non-matching elements removed

    Examples
    --------
    ::

        >>> import snakebids

    Filtering to get all ``subject='01'`` scans::

        >>> snakebids.filter_list(
        ...     {
        ...         'dir': ['AP','PA','AP','PA', 'AP','PA','AP','PA'],
        ...         'acq': ['98','98','98','98','99','99','99','99'],
        ...         'subject': ['01','01','02','02','01','01','02','02' ]
        ...     },
        ...     {'subject': '01'}
        ... ) == {
        ...     'dir': ['AP', 'PA', 'AP', 'PA'],
        ...     'acq': ['98', '98', '99', '99'],
        ...     'subject': ['01', '01', '01', '01']
        ... }
        True

    Filtering to get all ``acq='98'`` scans::

        >>> snakebids.filter_list(
        ...     {
        ...         'dir': ['AP','PA','AP','PA', 'AP','PA','AP','PA'],
        ...         'acq': ['98','98','98','98','99','99','99','99'],
        ...         'subject': ['01','01','02','02','01','01','02','02' ]
        ...     },
        ...     {'acq': '98'}
        ... ) == {
        ...     'dir': ['AP', 'PA', 'AP', 'PA'],
        ...     'acq': ['98', '98', '98', '98'],
        ...     'subject': ['01', '01', '02', '02']
        ... }
        True

    Filtering to get all ``dir=='AP'`` scans::

        >>> snakebids.filter_list(
        ...     {
        ...         'dir': ['AP','PA','AP','PA', 'AP','PA','AP','PA'],
        ...         'acq': ['98','98','98','98','99','99','99','99'],
        ...         'subject': ['01','01','02','02','01','01','02','02' ]
        ...     },
        ...     {'dir': 'AP'}
        ... ) == {
        ...     'dir': ['AP', 'AP', 'AP', 'AP'],
        ...     'acq': ['98', '98', '99', '99'],
        ...     'subject': ['01', '02', '01', '02']
        ... }
        True

    Filtering to get all ``subject='03'`` scans (i.e. no matches)::

        >>> snakebids.filter_list(
        ...     {
        ...         'dir': ['AP','PA','AP','PA', 'AP','PA','AP','PA'],
        ...         'acq': ['98','98','98','98','99','99','99','99'],
        ...         'subject': ['01','01','02','02','01','01','02','02' ]
        ...     },
        ...     {'subject': '03'}
        ... ) == {
        ...     'dir': [],
        ...     'acq': [],
        ...     'subject': []
        ... }
        True
    """
    keep_indices = []
    for key, val in wildcards.items():
        # get indices where wildcard exists
        if key not in zip_list.keys():
            continue
        indices = [i for i, v in enumerate(zip_list[key]) if v == val]
        if len(keep_indices) == 0:
            keep_indices = indices
        else:
            keep_indices = [x for x in indices if x in keep_indices]
    # Now we have the indices, so filter the lists
    if return_indices_only:
        return keep_indices
    return {
        key: [zip_list[key][i] for i in keep_indices] for key, val in zip_list.items()
    }


def get_filtered_ziplist_index(zip_list, wildcards, subj_wildcards):
    """Use this function when you have wildcards for a single scan instance,
    and want to know the index of that scan, amongst that subject's scan
    instances.

    Parameters
    ----------
    zip_list : dict
        lists for scans in a dataset, zipped to get each instance
    wildcards: dict
        wildcards for the single instance for querying it's index
    subj_wildcards: dict
        keys of this dictionary are used to pick out the subject/(session)
        from the wildcards

    Examples
    --------

    >>> import snakebids

    In this example, we have a dataset where with scans from two subjects,
    where each subject has ``dir-AP`` and ``dir-PA`` scans, along with
    ``acq-98`` and ``acq-99``:

    * ``sub-01_acq-98_dir-AP_dwi.nii.gz``
    * ``sub-01_acq-98_dir-PA_dwi.nii.gz``
    * ``sub-01_acq-99_dir-AP_dwi.nii.gz``
    * ``sub-01_acq-99_dir-PA_dwi.nii.gz``
    * ``sub-02_acq-98_dir-AP_dwi.nii.gz``
    * ``sub-02_acq-98_dir-PA_dwi.nii.gz``
    * ``sub-02_acq-99_dir-AP_dwi.nii.gz``
    * ``sub-02_acq-99_dir-PA_dwi.nii.gz``

    The ``zip_list`` produced by ``generate_inputs()`` is the set of entities
    that when zipped together, e.g. with ``expand(path, zip, **zip_list)``,
    produces the entity combinations that refer to each scan::

        {
            'dir': ['AP','PA','AP','PA', 'AP','PA','AP','PA'],
            'acq': ['98','98','98','98','99','99','99','99'],
            'subject': ['01','01','02','02','01','01','02','02']
        }

    The ``filter_list()`` function produces a subset of the entity
    combinations as a filtered zip list. This is used e.g. to get all the
    scans for a single subject.

    This ``get_filtered_ziplist_index()`` function performs ``filter_list()``
    twice:

    1. Using the ``subj_wildcards`` (e.g.: ``'subject': '{subject}'``) to get
       a subject/session-specific ``zip_list``.
    2. To return the indices from that list of the matching wildcards.

    In this example, if the wildcards parameter was::

        {'dir': 'PA', 'acq': '99', 'subject': '01'}

    Then the first (subject/session-specific) filtered list provides this zip
    list::

        {
            'dir': ['AP','PA','AP','PA'],
            'acq': ['98','98','99','99'],
            'subject': ['01','01','01','01']
        }

    which has 4 combinations, and thus are indexed from 0 to 3.

    The returned value would then be the index (or indices) that matches the
    wildcards. In this case, since the wildcards were
    ``{'dir': 'PA', 'acq': '99', 'subject':'01'}``, the return index is 3. ::

        >>> snakebids.get_filtered_ziplist_index(
        ...     {
        ...         'dir': ['AP','PA','AP','PA', 'AP','PA','AP','PA'],
        ...         'acq': ['98','98','98','98','99','99','99','99'],
        ...         'subject': ['01','01','02','02','01','01','02','02' ]
        ...     },
        ...     {'dir': 'PA', 'acq': '99', 'subject': '01'},
        ...     {'subject': '{subject}' }
        ... )
        3
    """
    # get the subject/(session) dict:
    subj_dict = {key: wildcards[key] for key in subj_wildcards.keys()}

    # now filter the list based on subj_wildcards
    zip_list_filtered = filter_list(zip_list, subj_dict)

    # get the index of the wildcard from this filtered list
    indices = filter_list(zip_list_filtered, wildcards, return_indices_only=True)
    if len(indices) == 1:
        return indices[0]
    return indices


def _read_bids_tags(bids_json=None):
    """Read the bids tags we are aware of from a JSON file. This is used
    specifically for compatibility with pybids, since some tag keys are
    different from how they appear in the file name, e.g. ``subject`` for
    ``sub``, and ``acquisition`` for ``acq``.

    Parameters
    ----------
    bids_json : str, optional
        Path to JSON file to use, if not specified will use
        ``bids_tags.json`` in the snakebids module.

    Returns
    -------
    dict:
        Dictionary of bids tags"""
    if bids_json is None:
        bids_json = Path(__file__).parent.resolve() / "bids_tags.json"
    with bids_json.open("r") as infile:
        bids_tags = json.load(infile)
    return bids_tags


def _generate_search_terms(participant_label=None, exclude_participant_label=None):
    search_terms = {}

    if participant_label is not None and exclude_participant_label is not None:
        raise ValueError(
            "Cannot define both participant_label and "
            "exclude_participant_label at the same time"
        )

    # add participant_label or exclude_participant_label to search terms (if
    # defined)
    # we make the subject key in search_terms a list so we can have both
    # participant_label and exclude_participant_label defined
    if participant_label is not None:
        if isinstance(participant_label, list):
            search_terms["subject"] = participant_label
        else:
            search_terms["subject"] = [participant_label]

    if exclude_participant_label is not None:
        # if multiple subjects to exclude, combine with with subj1|subj2|...
        if isinstance(exclude_participant_label, list):
            exclude_string = "|".join(exclude_participant_label)
        # if not, then string is the label itself
        else:
            exclude_string = exclude_participant_label
        search_terms["regex_search"] = True
        # regex to exclude subjects
        search_terms["subject"] = [f"^((?!({exclude_string})).)*$"]
    return search_terms


def generate_inputs(
    bids_dir,
    pybids_inputs,
    derivatives=False,
    search_terms=None,
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

    search_terms = _generate_search_terms(participant_label, exclude_participant_label)

    if os.path.exists(bids_dir):
        # generate inputs based on config
        layout = BIDSLayout(
            bids_dir,
            derivatives=derivatives,
            validate=False,
            indexer=BIDSLayoutIndexer(validate=False, index_metadata=False),
        )
    else:
        _logger.info(
            "bids_dir does not exist, skipping PyBIDS and using "
            "custom file paths only"
        )
        layout = None

    # this will populate input_path, input_lists, input_zip_lists, and
    # input_wildcards
    inputs_config_dict = _get_lists_from_bids(
        bids_layout=layout,
        pybids_inputs=pybids_inputs,
        limit_to=limit_to,
        **search_terms,
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
        inputs_config_dict["subjects"] = layout.get_subjects(**search_terms)
        inputs_config_dict["sessions"] = layout.get_sessions(**search_terms)

    if len(inputs_config_dict["sessions"]) == 0:
        inputs_config_dict["subj_wildcards"] = {"subject": "{subject}"}
    else:
        inputs_config_dict["subj_wildcards"] = {
            "subject": "{subject}",
            "session": "{session}",
        }

    return inputs_config_dict


def _process_path_override(input_path):
    """Glob wildcards from a path override and arrange into a zip list of
    matches, list of matches, and Snakemake wildcard for each wildcard.
    """
    wildcards = glob_wildcards(input_path)
    wildcard_names = list(wildcards._fields)

    if len(wildcard_names) == 0:
        _logger.warning("No wildcards defined in %s", input_path)

    input_wildcards = {}
    input_zip_lists = {}
    input_lists = {}

    for i, wildcard in enumerate(wildcard_names):
        input_zip_lists[wildcard] = wildcards[i]
        input_lists[wildcard] = list(set(wildcards[i]))
        input_wildcards[wildcard] = f"{{{wildcard}}}"
        if len(wildcards[i]) == 0:
            _logger.error("No matching files for %s", input_path)

    return input_zip_lists, input_lists, input_wildcards


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
    bids_tags = _read_bids_tags()
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
            ) = _process_path_override(input_path)
        else:
            paths = set()
            for img in bids_layout.get(
                **pybids_inputs[input_name]["filters"], **filters
            ):
                input_path = img.path
                for wildcard_name in pybids_inputs[input_name]["wildcards"]:

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
