import collections
from pathlib import Path, PosixPath, WindowsPath
from typing import Any, OrderedDict

from ruamel.yaml import YAML, Dumper


def get_yaml_io():
    """Return yaml loader/dumper configured for snakebids."""
    yaml = YAML()

    # Represent any PathLikes as str.
    def path2str(dumper: Dumper, data: Path):
        return dumper.represent_scalar(
            "tag:yaml.org,2002:str",
            str(data),
        )  # type: ignore

    def to_dict(dumper: Dumper, data: OrderedDict[Any, Any]):
        return dumper.represent_dict(dict(data))

    yaml.representer.add_representer(PosixPath, path2str)
    yaml.representer.add_representer(WindowsPath, path2str)
    yaml.representer.add_representer(collections.OrderedDict, to_dict)
    yaml.representer.add_representer(OrderedDict, to_dict)
    return yaml
