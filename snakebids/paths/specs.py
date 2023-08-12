import importlib_resources as impr
import more_itertools as itx
import yaml
from typing_extensions import NotRequired, TypeAlias, TypedDict

from snakebids.paths import resources


class BidsPathEntitySpec(TypedDict):
    entity: str
    tag: NotRequired[str]
    dir: NotRequired[bool]


BidsPathSpec: TypeAlias = "list[BidsPathEntitySpec]"


def _find_entity(spec: BidsPathSpec, entity: str):
    return itx.one(item for item in spec if item["entity"] == entity)


def v0_0_0(subject_dir: bool = True, session_dir: bool = True) -> BidsPathSpec:
    """Get the v0.0.0 BidsPathSpec

    This is the legacy spec used since the beginning of snakebids.

    Formatted as ``sub-*/ses-*/{datatype}/sub-*_ses-*_task-*_acq-*_ce-*_rec-*_dir-*_run\
-*_mod-*_echo-*_hemi-*_space-*_res-*_den-*_label-*_desc-*_..._{suffix}{.ext}``

    Parameters
    ----------
    subject_dir
        If False, downstream path generator will not include the subject dir
        `sub-{subject}/*`
    session_dir : bool, optional
        If False, downstream path generator will not include the session dir
        `*/ses-{session}/*`
    """
    spec = yaml.safe_load(impr.files(resources).joinpath("spec.0.0.0.yaml").read_text())
    if not subject_dir:
        _find_entity(spec, "subject")["dir"] = False

    if not session_dir:
        _find_entity(spec, "session")["dir"] = False

    return spec
