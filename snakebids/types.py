from __future__ import annotations

from typing_extensions import TypeAlias, TypedDict


class InputConfig(TypedDict, total=False):
    """Configuration passed in snakebids.yaml file"""

    filters: dict[str, str | bool | list[str]]
    wildcards: list[str]
    custom_path: str


InputsConfig: TypeAlias = "dict[str, InputConfig]"
