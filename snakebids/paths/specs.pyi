# This stub file is automatically generated
# It can be updated using::
#
#      poetry run poe update_bids

from ._utils import BidsPathSpec

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

def v0_11_0(subject_dir: bool = True, session_dir: bool = True) -> BidsPathSpec:
    """Spec corresponding to `BIDS v1.9.0`_.

    Significantly expanded from the v0.0.0 spec, now including long names for every
    relevant entity. In addition to the official spec, it includes `from` and `to`
    entities intended for transformations. Unknown entities are placed just before desc,
    so that the description entity is always last.

    .. _BIDS v1.9.0: https://bids-specification.readthedocs.io/en/v1.9.0/

    Formatted as::

        sub-{subject}/ses-{session}/{datatype}/{prefix}_sub-{subject}_ses-{session}_
        sample-{sample}_task-{task}_tracksys-{tracksys}_acq-{acquisition}_
        ce-{ceagent}_stain-{staining}_trc-{tracer}_rec-{reconstruction}_
        dir-{direction}_run-{run}_mod-{modality}_echo-{echo}_flip-{flip}_
        inv-{inversion}_mt-{mt}_proc-{processed}_part-{part}_space-{space}_
        atlas-{atlas}_seg-{segmentation}_hemi-{hemisphere}_res-{resolution}_
        den-{density}_roi-{roi}_from-{from}_to-{to}_split-{split}_
        recording-{recording}_chunk-{chunk}_model-{model}_subset-{subset}_
        label-{label}_..._desc-{description}_{suffix}{extension}


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
    """Spec corresponding to `BIDS v1.9.0`_.

    Significantly expanded from the v0.0.0 spec, now including long names for every
    relevant entity. In addition to the official spec, it includes `from` and `to`
    entities intended for transformations. Unknown entities are placed just before desc,
    so that the description entity is always last.

    .. _BIDS v1.9.0: https://bids-specification.readthedocs.io/en/v1.9.0/

    Formatted as::

        sub-{subject}/ses-{session}/{datatype}/{prefix}_sub-{subject}_ses-{session}_
        sample-{sample}_task-{task}_tracksys-{tracksys}_acq-{acquisition}_
        ce-{ceagent}_stain-{staining}_trc-{tracer}_rec-{reconstruction}_
        dir-{direction}_run-{run}_mod-{modality}_echo-{echo}_flip-{flip}_
        inv-{inversion}_mt-{mt}_proc-{processed}_part-{part}_space-{space}_
        atlas-{atlas}_seg-{segmentation}_hemi-{hemisphere}_res-{resolution}_
        den-{density}_roi-{roi}_from-{from}_to-{to}_split-{split}_
        recording-{recording}_chunk-{chunk}_model-{model}_subset-{subset}_
        label-{label}_..._desc-{description}_{suffix}{extension}


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
