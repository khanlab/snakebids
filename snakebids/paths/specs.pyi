# This stub file is automatically generated
# It can be updated using::
#
#      poetry run poe update_bids

from .utils import BidsPathSpec

LATEST: str

def v0_0_0(subject_dir: bool = True, session_dir: bool = True) -> BidsPathSpec:
    """Get the v0.0.0 BidsPathSpec.

    This spec alone equips :func:`~snakebids.bids` with 2 extra arguments:
    ``include_subject_dir`` and ``include_session_dir``. These default to ``True``, but
    if set ``False``, remove the subject and session dirs respectively from the output
    path. For future specs, this behaviour should be achieved by modifying the spec and
    generating a new :func:`~snakebids.bids` function

    Formatted as::

        sub-{subject}/ses-{session}/{datatype}/{prefix}_sub-{subject}_ses-{session}_
        task-{task}_acq-{acq}_ce-{ce}_rec-{rec}_dir-{dir}_run-{run}_mod-{mod}_
        echo-{echo}_hemi-{hemi}_space-{space}_res-{res}_den-{den}_label-{label}_
        desc-{desc}_..._{suffix}{extension}


    Parameters
    ----------
    subject_dir
        If False, downstream path generator will not include the subject dir
        `sub-{subject}/*`
    session_dir : bool, optional
        If False, downstream path generator will not include the session dir
        `*/ses-{session}/*`
    """
    ...

def latest(subject_dir: bool = True, session_dir: bool = True) -> BidsPathSpec:
    """Get the v0.0.0 BidsPathSpec.

    This spec alone equips :func:`~snakebids.bids` with 2 extra arguments:
    ``include_subject_dir`` and ``include_session_dir``. These default to ``True``, but
    if set ``False``, remove the subject and session dirs respectively from the output
    path. For future specs, this behaviour should be achieved by modifying the spec and
    generating a new :func:`~snakebids.bids` function

    Formatted as::

        sub-{subject}/ses-{session}/{datatype}/{prefix}_sub-{subject}_ses-{session}_
        task-{task}_acq-{acq}_ce-{ce}_rec-{rec}_dir-{dir}_run-{run}_mod-{mod}_
        echo-{echo}_hemi-{hemi}_space-{space}_res-{res}_den-{den}_label-{label}_
        desc-{desc}_..._{suffix}{extension}


    Parameters
    ----------
    subject_dir
        If False, downstream path generator will not include the subject dir
        `sub-{subject}/*`
    session_dir : bool, optional
        If False, downstream path generator will not include the session dir
        `*/ses-{session}/*`
    """
    ...
