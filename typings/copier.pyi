from os import PathLike
from typing import Any, AnyStr, Literal, Mapping, TypedDict

from git import Sequence
from typing_extensions import TypeAlias, Unpack

AnyPath: TypeAlias = "AnyStr | PathLike[AnyStr]"

class CopierArgs(TypedDict, total=False):
    answers_file: AnyPath[str] | AnyPath[bytes]
    vcs_ref: str | None
    exclude: Sequence[str]
    use_prereleases: bool
    skip_if_exists: Sequence[str]
    cleanup_on_error: bool
    defaults: bool
    user_defaults: dict[str, Any]
    overwrite: bool
    pretend: bool
    quiet: bool
    conflict: Literal["inline", "rej"]
    context_lines: int
    unsafe: bool
    skip_answered: bool

def run_copy(
    src_path: str,
    dst_path: AnyPath[str] | AnyPath[bytes] = ...,
    data: Mapping[str, Any] | None = ...,
    **kwargs: Unpack[CopierArgs],
) -> None: ...
