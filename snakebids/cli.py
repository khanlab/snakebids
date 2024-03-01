from __future__ import annotations

import argparse
import warnings
from typing import Any, Mapping

from snakebids.types import InputsConfig


def add_dynamic_args(
    parser: argparse.ArgumentParser,
    parse_args: Mapping[str, Any],
    pybids_inputs: InputsConfig,
) -> None:
    """Do nothing.

    Originally added --filter-<comp> and --wildcards-<comp> argumets to the CLI. Kept
    as a placeholder for apps that relied on it for generating documentation. This
    functionality is now native to `SnakeBidsApp`.
    """
    warnings.warn(
        "add_dynamic_args() is deprecated and no longer has any effect. Its function "
        "is now provided natively by `SnakeBidsApp`. It will be removed in an upcoming "
        "release",
        stacklevel=2,
    )
