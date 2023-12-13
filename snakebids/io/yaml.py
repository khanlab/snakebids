from pathlib import Path, PosixPath, WindowsPath

from ruamel.yaml import YAML, Dumper


def get_yaml_io():
    """Return yaml loader/dumper configured for snakebids."""
    yaml = YAML()

    # Represent any PathLikes as str.
    def path2str(dumper: Dumper, data: Path):
        return dumper.represent_scalar(
            "tag:yaml.org,2002:str",
            str(data),
        )

    yaml.representer.add_representer(PosixPath, path2str)
    yaml.representer.add_representer(WindowsPath, path2str)
    return yaml
