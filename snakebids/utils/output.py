"""Tools to manage generation and interconversion of bidsapps and snakemake outputs."""

import hashlib
import itertools as it
import json
import os
import shutil
import time
from collections import OrderedDict
from pathlib import Path, PosixPath, WindowsPath
from typing import Iterable, List, Union

import yaml
from colorama import Fore
from progress.bar import IncrementalBar
from typing_extensions import Literal

from snakebids.exceptions import RunError

Mode = Union[Literal["workflow"], Literal["bidsapp"]]


def prepare_output(
    src: Path, outputdir: Path, mode: Mode, force_conversion: bool = False
):
    """Ensure output directory is in the correct mode and is ready for snakemake to run

    Checks for existing output at the directory, converting to the appropriate mode if
    necessary. If running in workflow mode, the src snakemake directory is copied over
    to the output, and snakemake output will be put in results. In bidsapp mode, the
    snakemake results will be put in the top level of the output directory without any
    extra snakemake files. Creates a .snakebids file to track output mode and other
    workflow information. Raises exceptions when attempting to convert from workflow
    mode to bidsapp mode without force, when the output directory already has contents
    but no .snakebids file to identify it, or when attempting to convert snakemake
    results that, themselves, have a results folder

    Parameters
    ----------
    src : Path
        Path to snakebids app
    outputdir : Path
        Path to output
    mode : "bidsapp" or "workflow"
        Mode in which to run
    force_conversion : bool
        Force conversion from workflow to bidsapp
        mode. Defaults to False.

    Returns
    -------
    Path
        Path to new root folder (output for bidsapp, output/results for workflow)

    Raises
    ------
    RunError
        Raised when attempting to convert from workflow mode to bidsapp mode
        without force, when the output directory already has contents but no
        .snakebids file to identify it, or when attempting to convert snakemake
        results that, themselves, have a results folder

    """
    if mode not in ["workflow", "bidsapp"]:
        raise RunError(
            f"Requested unsupported output mode: {mode}.\n"
            'Please select between "workflow" and "bidsapp"'
        )

    # Look for .snakebids file. If the outputdir doesn't yet exist, we'll get None.
    # If it does exist but there's no .snakebids file, an error will be raised.
    snakebids_file = _get_snakebids_file(outputdir)
    outputdir.mkdir(exist_ok=True)

    if not snakebids_file:
        write_output_mode(outputdir / ".snakebids", mode)

        if mode == "bidsapp":
            return outputdir

        return _copy_snakemake_app(src, outputdir)

    with snakebids_file.open("r") as f:
        snakebids_data = json.load(f)
    root = _convert_output(
        snakebids_data["mode"], mode, src, outputdir, force_conversion
    )
    write_output_mode(outputdir / ".snakebids", mode)
    return root


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


def retrofit_output(output: Path, config_files: Iterable[Path]):
    """Convert legacy snakebids output to bidsapp mode.

    Expects a directory containing previous snakebids outputs. This should contain one
    or more config files, as specified in the config_files parameter. If a config
    directory resides in the output, it should only contain specified config files. If
    extra files are found, an error will be thrown. All config files will be deleted,
    and a new .snakebids file will be created.

    Parameters
    ----------
    output : Path
        Path of output directory to be modified
    config_files : Iterable[Path]
        List of paths of config files.

    Returns
    -------
    boolean
        Returns True if successful, False if User declined operation in interactive
        prompt

    Raises
    ------
    RunError
        Raised if existing .snakebids file found, if unrecognized files are found in the
        config folder, or if no config files are provided.
    """
    config_files = [*config_files]
    if (output / ".snakebids").exists():
        raise RunError(f".snakebids file already found at {output}")
    if not config_files:
        raise RunError(f"No config files found in {output}. Cannot perform retrofit.")
    if (output / "config").exists() and (output / "config").is_dir():
        to_delete = it.chain(config_files, [output / "config"])
        unknown_files = {*(output / "config").iterdir()} - {*config_files}
        if len(unknown_files) > 0:
            raise RunError(
                f"Unrecognized files found in config folder ({output/'config'}",
                _format_path_list(unknown_files),
            )
    else:
        to_delete = config_files

    print(f"Converting {output} into bidsapp format.\n")

    if not _remove_all(to_delete, confirm=True):
        return False
    write_output_mode(output / ".snakebids", "bidsapp")
    return True


def _convert_output(start: Mode, end: Mode, src: Path, output: Path, force=False):
    """Convert existing output between bidsapp and workflow mode

    Does nothing if start and end are the same. Because of the potential loss of
    information (all the workflow files), this will only convert from bidsapp to
    workflow mode if force is True, otherwise it raises an exception. Returns the
    root folder of the converted dataset (ouput if bidsapp mode, results if workflow)

    Parameters
    ----------
    start : "bidsapp" or "workflow"
        Current format of the output
    end : "bidsapp" or "workflow"
        Desired format of the output
    src : Path
        Path to the snakebids app being run
    output : Path
        Output to transform
    force : bool
        Force conversions from workflow to bidsapp mode.
        Defaults to False.

    Returns
    -------
    Path
        Path to the root folder (output/results if workflow, output if bidsapp)

    Raises
    ------
    RunError
        Raised if conversion from workflow to bidsapp attempted without force

    """
    if start == end:
        if end == "workflow":
            return output / "results"
        return output

    # Convert to workflow mode
    if end == "workflow":
        print(
            f"{Fore.GREEN}Found bidsapp output at {Fore.RESET}{output.resolve()}\n\n"
            f"{Fore.GREEN}Converting to workflow mode...{Fore.RESET}\n"
        )
        results = _check_for_results_folder(output)
        results.mkdir()
        for f in output.iterdir():
            if f != results:
                shutil.move(f, results / f.name)

        # Move .snakebids file back to the top level
        shutil.move(results / ".snakebids", output / ".snakebids")
        _copy_snakemake_app(src, output, False)
        return results

    # Convert to Bidsapp mode
    if not force:
        raise RunError(
            f"You are attempting to convert a preexisting output ({output}) "
            "from workflow mode to bidsapp mode. This will result in the loss "
            "of all workflow files and custom configs. If you are sure you "
            "wish to do this, run this command again with --force-conversion. "
            "Otherwise, run the command with --workflow-mode to maintain the current "
            "output mode."
        )
    # Check if results folder/file exists within the results folder. We don't
    # need its output, just its exception.
    _check_for_results_folder(output / "results")

    # Delete everything in the output folder except for .snakebids and results
    _remove_all(
        f
        for f in output.iterdir()
        if f not in [output / ".snakebids", output / "results", output / ".snakemake"]
    )

    for f in (output / "results").iterdir():
        shutil.move(f, output / f.name)

    (output / "results").rmdir()
    return output


def _remove_all(paths: Iterable[Path], confirm: bool = False):
    if confirm:
        paths = [*paths]
        print(
            f"\t{Fore.YELLOW}The following files and folders will be DELETED:\n"
            f"{Fore.RESET}",
            _format_path_list(paths),
        )
        user_response = input("Would you like to continue? [yes,NO]")
        if user_response.lower() != "yes":
            return False

    dirs: List[Path] = []
    for path in paths:
        if path.is_dir():
            dirs.append(path)
        if path.is_symlink():
            path.unlink()
        if path.is_file():
            os.remove(path)
    for d in dirs:
        shutil.rmtree(d)
    return True


def _format_path_list(paths: Iterable[Path]):
    return "\t\t- " + "\n\t\t- ".join(str(p) for p in paths)


def _copy_snakemake_app(src: Path, dest: Path, create_results: bool = True):
    """Copies snakemake app from src to dest, skipping config and results directories

    Creates an empty results and config folder, and returns the path to the results
    folder.

    Parameters
    ----------
    src : Path
        Directory to copy from
    dest : Path
        Directory to copy to
    create_results: bool
        If True, create an empty results folder (Default value = True)

    Returns
    -------
    Path
        Path of results folder in new snakemake location

    """
    # shutils.copytree makes it hard to exclude just root level folders, so we do this
    # manually
    file_list = [
        *it.chain(
            # All root level items
            f.relative_to(src)
            for f in Path(src).iterdir()
            if f
            not in [
                src / ".snakebids",
                src / "config",
                src / "results",
                src / ".snakemake",
            ]
        )
    ]

    # Iterate through files
    progress_bar = IncrementalBar("Copying snakebids app", max=len(file_list))
    for file in file_list:
        if (src / file).is_dir():
            shutil.copytree(src / file, dest / file)
        # Copy over files, making parents as necessary
        elif (src / file).is_file():
            (dest / file).parent.mkdir(parents=True, exist_ok=True)
            shutil.copy(src / file, dest / file, follow_symlinks=False)
        progress_bar.next()

    if create_results:
        (dest / "results").mkdir()
    (dest / "config").mkdir()
    return dest / "results"


def _get_snakebids_file(outputdir: Path):
    """Ensure populated dir contains .snakebids file, retrieving it if it does.

    First checks if outputdir doesn't exist or is completely empty, returing None if so.
    If it does have data, it checks for a .snakebids file, returning its Path if found.
    If no .snakebids file is found, it raises an exception.

    Parameters
    ----------
    outputdir : Path
        Directory to check.

    Returns
    -------
    Path or None
        None if output dir is nonexistant or empty, otherwise the path
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
    if len([*outputdir.iterdir()]) == 0:
        return None

    # If it's not empty, is there a .snakebids file?
    if (outputdir / ".snakebids").exists():
        return outputdir / ".snakebids"

    # We have an occupied directory without a .snakebids file, so we have no idea
    # what's there.
    raise RunError(
        f"Output dir `{outputdir.resolve()}` exists, but `.snakebids` file "
        "not found. Please specify either a new directory, or a ",
        "directory where you've previously run this Snakebids app.",
    )


def _check_for_results_folder(root: Path):
    """Check folder for results folder

    Raises exception if it does exist, otherwise returns the name of the folder it
    looked for.

    Parameters
    ----------
    root : Path
        Folder in which to search

    Returns
    -------
    Path
        Name of folder searched for.

    Raises
    ------
    RunError
        Raised if results folder found

    """
    results = root / "results"
    if results.exists():
        raise RunError(
            "Cannot convert output format as a results folder or file already "
            f"exists in your data directory ({results}). Please rename or remove "
            "this item."
        )
    return results


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
