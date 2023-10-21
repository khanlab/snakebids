from __future__ import annotations

from typing import List

import importlib_resources as impr
import more_itertools as itx
from typing_extensions import NotRequired, TypeAlias, TypedDict

from snakebids.io.yaml import get_yaml_io
from snakebids.paths import resources


class BidsPathEntitySpec(TypedDict):
    """Interface for BIDS path specification."""

    entity: str
    """Entity full name"""

    tag: NotRequired[str]
    """Short entity name, as appears in the path"""

    dir: NotRequired[bool]
    """If true, a directory with the entity-value pair is created"""


def _find_entity(spec: BidsPathSpec, entity: str):
    return itx.one(item for item in spec if item["entity"] == entity)


BidsPathSpec: TypeAlias = List[BidsPathEntitySpec]
"""List of :class:`BidsPathEntitySpec`, defining the order of entities in a bids path"""


def v0_0_0(subject_dir: bool = True, session_dir: bool = True) -> BidsPathSpec:
    r"""Get the v0.0.0 BidsPathSpec.

    This spec alone equips :func:`~snakebids.bids` with 2 extra arguments:
    ``include_subject_dir`` and ``include_session_dir``. These default to ``True``, but
    if set ``False``, remove the subject and session dirs respectively from the output
    path. For future specs, this behaviour should be achieved by modifying the spec and
    generating a new :func:`~snakebids.bids` function

    Formatted as::

        sub-{sub}/ses-{ses}/{datatype}/\
            sub-{sub}_ses-{ses}_task-{task}_acq-{acq}_\
            ce-{ce}_rec-{rec}_dir-{dir}_run-{run}_mod-{mod}_\
            echo-{echo}_hemi-{hemi}_space-{space}_res-{res}_\
            den-{den}_label-{label}_desc-{desc}_..._{suffix}{.ext}

    Parameters
    ----------
    subject_dir
        If False, downstream path generator will not include the subject dir
        `sub-{subject}/*`
    session_dir : bool, optional
        If False, downstream path generator will not include the session dir
        `*/ses-{session}/*`
    """
    spec = get_yaml_io().load(
        impr.files(resources).joinpath("spec.0.0.0.yaml").read_text()
    )
    if not subject_dir:
        _find_entity(spec, "subject")["dir"] = False

    if not session_dir:
        _find_entity(spec, "session")["dir"] = False

    return spec
