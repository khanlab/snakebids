from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING, TypedDict

from typing_extensions import TypeAlias

if TYPE_CHECKING:
    from argparse import _ArgumentGroup

__all__ = ["ArgumentGroups", "SnakebidsConfig"]

ArgumentGroups: TypeAlias = "dict[str, _ArgumentGroup]"


class SnakebidsConfig(TypedDict):
    """Interface of default arguments parsed by {func}`snakebids.bidsapp.parser`."""

    bids_dir: Path
    output_dir: Path
    analysis_level: str
    participant_label: list[str] | None
    exclude_participant_label: list[str] | None
    derivatives: list[str] | None
