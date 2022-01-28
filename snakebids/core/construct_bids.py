"""Utilities for converting Snakemake apps to BIDS apps."""

from collections import OrderedDict
from pathlib import Path
from typing import Optional, Union


# pylint: disable=too-many-arguments
def bids(
    root: Union[str, Path] = None,
    datatype: str = None,
    prefix: str = None,
    suffix: str = None,
    subject: str = None,
    session: str = None,
    include_subject_dir: bool = True,
    include_session_dir: bool = True,
    **entities: str,
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
    Path
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

    # replace underscores in keys (needed so that users can use reserved
    # keywords by appending a _)
    entities = {k.replace("_", ""): v for k, v in entities.items()}

    # strict ordering of bids entities is specified here:
    # pylint: disable=unsubscriptable-object
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

    # Form all entities for filename as a list, and join with "_". Any undefined
    # entities will be `None` and will be filtered out.
    filename: str = "_".join(
        filter(
            None,
            [
                # Put the prefix before anything else
                prefix,
                # Add in subject and session
                f"sub-{subject}" if subject else None,
                f"ses-{session}" if session else None,
                # Iterate through all other entities and add as "key-value"
                *(f"{key}-{val}" for key, val in order.items() if val is not None),
                # Put the suffix last
                suffix,
            ],
        )
    )

    # If all entities were `None`, the list will be empty and filename == ""
    if not filename:
        return ""

    # Form folder using list similar to filename, above. Filter out Nones, and convert
    # to Path.
    folder = Path(
        *filter(
            None,
            [
                str(root) if root else None,
                f"sub-{subject}" if subject and include_subject_dir else None,
                f"ses-{session}" if session and include_session_dir else None,
                datatype,
            ],
        )
    )

    return str(folder / filename)


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
