"""Tools to manage generation and interconversion of bidsapps and snakemake outputs."""

import hashlib
import json
import time
from collections import OrderedDict
from pathlib import Path, PosixPath, WindowsPath
from typing import Dict, Union

import more_itertools as itx
import yaml
from typing_extensions import Literal

from snakebids.exceptions import RunError

Mode = Union[Literal["workflow"], Literal["bidsapp"]]


def prepare_bidsapp_output(outputdir: Path, force_output: bool = False):
    """Ensure output directory is in the correct mode and is ready for snakemake to run

    Checks for existing output at the directory, creating it if necessary. Creates a
    .snakebids file to track output mode and other workflow information. If outputdir
    already has contents but no .snakebids file, it raises an exception unless
    force_output == True.

    Parameters
    ----------
    outputdir : Path
        Path to output
    force_output : bool
        Force output in a new output directory that already has contents. Defaults to
        False.

    Returns
    -------
    Path
        Path to new root folder (output for bidsapp, output/results for workflow)

    Raises
    ------
    RunError
        Raised when the output directory already has contents but no .snakebids file to
        identify it

    """
    # Look for .snakebids file. If the outputdir doesn't yet exist, we'll get None.
    # If it does exist but there's no .snakebids file, an error will be raised.
    try:
        _get_snakebids_file(outputdir)
    except RunError as err:
        if not force_output:
            raise err

    outputdir.mkdir(exist_ok=True)

    write_output_mode(outputdir / ".snakebids", "bidsapp")


def write_output_mode(dotfile: Path, mode: Mode):
    """Write output mode to .snakebids

    Parameters
    ----------
    dotfile: Path
        Path to .snakebids file to be written
    mode: Mode
        Mode to write: either "bidsapp" or "workflow"
    """
    if dotfile.exists():
        with dotfile.open("r") as f:
            data = json.load(f)
    else:
        data = {}
    data["mode"] = mode
    with (dotfile).open("w") as f:
        json.dump(data, f)


def _get_snakebids_file(outputdir: Path):
    """Ensure populated dir contains .snakebids file, retrieving it if it does.

    First checks if outputdir doesn't exist or is completely empty, returing None if so.
    If it does have data, it checks for a .snakebids file, returning its contents if
    found. If no .snakebids file is found, it raises an exception.

    Parameters
    ----------
    outputdir : Path
        Directory to check.

    Returns
    -------
    Dict or None
        None if output dir is nonexistant or empty, otherwise the contents
        of the .snakebids file

    Raises
    ------
    RunError
        Raised if outputdir contains contents but no .snakebids file

    """
    # Check if outputdir exits
    if not outputdir.exists():
        return None

    # If it does exist, is it empty?
    if itx.ilen(outputdir.iterdir()) == 0:
        return None

    # If it's not empty, is there a .snakebids file?
    if (outputdir / ".snakebids").exists():
        malformed_err = RunError(
            f"Found malformed .snakebids file in `{outputdir.resolve()}. Please"
            "remove this file and check the integrity of previous outputs."
        )
        with (outputdir / ".snakebids").open("r") as f:
            try:
                snakebids_data: Dict[str, str] = json.load(f)
            except json.JSONDecodeError as err:
                raise malformed_err from err
        if "mode" not in snakebids_data:
            raise malformed_err

        return snakebids_data

    # We have an occupied directory without a .snakebids file, so we have no idea
    # what's there.
    raise RunError(
        f"Output dir `{outputdir.resolve()}` exists, but `.snakebids` file "
        "not found. Please specify either a new directory, or a ",
        "directory where you've previously run this Snakebids app.",
    )


def get_time_hash():
    """currently unused"""

    time_hash = hashlib.sha1()
    time_hash.update(str(time.time()).encode("utf-8"))
    return time_hash.hexdigest()[:8]


def write_config_file(config_file: Path, data: dict, force_overwrite: bool = False):
    if (config_file.exists()) and not force_overwrite:
        raise RunError(
            f"A config file named {config_file.name} already exists:\n"
            f"\t- {config_file.resolve()}\n"
            "Please move or rename either the existing or incoming config."
        )
    config_file.parent.mkdir(exist_ok=True)

    # TODO: copy to a time-hashed file for provenance purposes?
    #       unused as of now..
    # time_hash = get_time_hash()

    with open(config_file, "w", encoding="utf-8") as f:
        # write either as JSON or YAML
        if config_file.suffix == ".json":
            json.dump(data, f, indent=4)
            return

        # if not json, then should be yaml or yml

        # this is needed to make the output yaml clean
        yaml.add_representer(
            OrderedDict,
            lambda dumper, data: dumper.represent_mapping(
                "tag:yaml.org,2002:map", data.items()
            ),
        )

        # Represent any PathLikes as str.
        def path2str(dumper, data):
            return dumper.represent_scalar("tag:yaml.org,2002:str", str(data))

        yaml.add_representer(PosixPath, path2str)
        yaml.add_representer(WindowsPath, path2str)

        yaml.dump(data, f, default_flow_style=False, sort_keys=False)
