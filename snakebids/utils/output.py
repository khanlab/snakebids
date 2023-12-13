"""Tools to manage generation and interconversion of bidsapps and snakemake outputs."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Literal

import more_itertools as itx

from snakebids.exceptions import RunError

Mode = Literal["workflow", "bidsapp"]


def prepare_bidsapp_output(outputdir: Path, force_output: bool = False) -> None:
    """Ensure output directory is in the correct mode and is ready for snakemake to run.

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
    except RunError:
        if not force_output:
            raise

    outputdir.mkdir(exist_ok=True)

    write_output_mode(outputdir / ".snakebids", "bidsapp")


def write_output_mode(dotfile: Path, mode: Mode) -> None:
    """Write output mode to .snakebids.

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


def _get_snakebids_file(outputdir: Path) -> dict[str, str] | None:
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
            f"Found malformed .snakebids file in `{outputdir.resolve()}. Please "
            "remove this file and check the integrity of previous outputs."
        )
        with (outputdir / ".snakebids").open("r") as f:
            try:
                snakebids_data: dict[str, str] = json.load(f)
            except json.JSONDecodeError as err:
                raise malformed_err from err
        if "mode" not in snakebids_data:
            raise malformed_err

        return snakebids_data

    # We have an occupied directory without a .snakebids file, so we have no idea
    # what's there.
    msg = (
        f"Output dir `{outputdir.resolve()}` exists, but it doesn't look like this "
        "app has been run there before. If you're sure you got the directory correct, "
        "run the app again using `--force-output`"
    )
    raise RunError(
        msg,
    )
