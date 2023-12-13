from __future__ import annotations

import errno
import json
from pathlib import Path
from typing import TYPE_CHECKING, Any

from snakebids.io.yaml import get_yaml_io

if TYPE_CHECKING:
    from _typeshed import StrPath


def write_config(
    config_file: StrPath, data: dict[str, Any], force_overwrite: bool = False
) -> None:
    """Write provided data as yaml or json to provided path.

    Output type is decided based on file suffix: .json -> JSON, .yaml,.yml -> YAML

    Parameters
    ----------
    config_file
        Path of yaml file
    data
        Data to format
    force_overwrite
        If True, force overwrite of already existing files, otherwise error out
    """
    config_file = Path(config_file)
    if (config_file.exists()) and not force_overwrite:
        err = FileExistsError(
            errno.EEXIST, f"'{config_file}' already exists", str(config_file)
        )
        raise err
    config_file.parent.mkdir(exist_ok=True)

    if config_file.suffix == ".json":
        with open(config_file, "w", encoding="utf-8") as f:
            json.dump(data, f, indent=4)
            return

    get_yaml_io().dump(data, config_file)
