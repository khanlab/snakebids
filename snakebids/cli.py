import argparse
import os
import pathlib
import re
from typing import Any, Dict, List

import attr
import snakemake

# We define Path here in addition to pathlib to put both variables in globals()
# This way, users specifying a path type in their config.yaml can indicate
# either Path or pathlib.Path
Path = pathlib.Path


# pylint: disable=too-few-public-methods,
class KeyValue(argparse.Action):
    """Class for accepting key=value pairs in argparse"""

    # Constructor calling
    def __call__(self, parser, namespace, values, option_string=None):
        setattr(namespace, self.dest, {})

        for value in values:
            # split it into key and value
            key, value = value.split("=")
            # assign into dictionary
            getattr(namespace, self.dest)[key] = value


# pylint: disable=too-few-public-methods
class SnakemakeHelpAction(argparse.Action):
    """Class for printing snakemake usage in argparse"""

    def __call__(self, parser, namespace, values, option_string=None):
        snakemake.main(["-h"])


# pylint: disable=missing-class-docstring
@attr.frozen
class SnakebidsArgs:
    """Arguments for Snakebids App

    Organizes the various options available for a generic Snakebids App, and store
    project specific arguments in a dict. Snakemake args are to be put in a list

    Attributes
    ----------
    workflow_mode : bool
        Indicates whether to run in workflow mode. False maps to bidsapp mode
    force : bool
        Force conversion from workflow apps to bidsapp apps, resulting in potential
        data loss
    outputdir : Path
        Directory to place outputs
    retrofit : Path
        Convert a legacy snakebids app to the new style
    snakemake_args : list of strings
        Arguments to pass on to Snakemake
    args_dict : Dict[str, Any]
        Contains all the snakebids specific args. Meant to contain custom user args
        defined in config, as well as dynamic --filter-xx and --wildcard-xx args.
        These will eventually be printed in the new config.
    """

    workflow_mode: bool
    force: bool
    outputdir: Path
    retrofit: bool
    snakemake_args: List[str]
    args_dict: Dict[str, Any]


def create_parser(include_snakemake=False):
    """Generate basic Snakebids Parser

    Includes the standard Snakebids arguments.
    """

    # The snakemake parser functionality did not seem to be implemented in the original
    # factoring. I left the logic here, but it should probably be removed if it's not
    # needed.
    if include_snakemake:
        # get snakemake parser
        smk_parser = snakemake.get_argument_parser()

        # create parser
        parser = argparse.ArgumentParser(
            description="Snakebids helps build BIDS Apps with Snakemake",
            add_help=False,
            parents=[smk_parser],
        )
    else:
        parser = argparse.ArgumentParser(
            description="Snakebids helps build BIDS Apps with Snakemake"
        )

    standard_group = parser.add_argument_group(
        "STANDARD", "Standard options for all snakebids apps"
    )

    standard_group.add_argument(
        "--workflow-mode",
        "--workflow_mode",
        "-W",
        action="store_true",
        help=(
            "Run Snakebids in Workflow mode. This will cause the entire Snakebids "
            "app, except for the results folder, to be copied into your "
            "output_dir. Snakemake will be run in this new child app, and results "
            "will be put in output_dir/results."
        ),
    )

    # We use -x as the alias because both -f and -F are taken by snakemake
    standard_group.add_argument(
        "--force-conversion",
        "--force_conversion",
        "-x",
        action="store_true",
        help=(
            "Force snakebids to convert a workflow output to a bidsapp output. "
            "conversion will delete all the files in the workflow snakebids app."
        ),
    )

    standard_group.add_argument(
        "--retrofit",
        action="store_true",
        help=(
            "Convert a legacy Snakebids output (Snakebids 0.3.x or lower) into "
            "bidsapp mode. This will delete any config files currently in the "
            "output."
        ),
    )

    # add option for printing out snakemake usage
    standard_group.add_argument(
        "--help-snakemake",
        "--help_snakemake",
        nargs=0,
        action=SnakemakeHelpAction,
        help=(
            "Options to Snakemake can also be passed directly at the "
            "command-line, use this to print Snakemake usage"
        ),
    )
    return parser


def add_dynamic_args(
    parser: argparse.ArgumentParser,
    parse_args: Dict[str, Dict[str, str]],
    pybids_inputs: Dict[str, Dict[str, Dict[str, Any]]],
):
    # create parser group for app options
    app_group = parser.add_argument_group("SNAKEBIDS", "Options for snakebids app")

    # update the parser with config options
    for name, arg in parse_args.items():
        # Convert type annotations from strings to class types
        # We first check that the type annotation is, in fact,
        # a str to allow the edge case where it's already
        # been converted
        if "type" in arg and isinstance(arg["type"], str):
            try:
                arg["type"] = globals()[arg["type"]]
            except KeyError as err:
                raise TypeError(
                    f"{arg['type']} is not available " + f"as a type for {name}"
                ) from err

        names = _make_underscore_dash_aliases(name)
        app_group.add_argument(*names, **arg)

    # general parser for
    # --filter_{input_type} {key1}={value1} {key2}={value2}...
    # create filter parsers, one for each input_type
    filter_opts = parser.add_argument_group(
        "BIDS FILTERS",
        "Filters to customize PyBIDS get() as key=value pairs",
    )

    for input_type in pybids_inputs.keys():
        argnames = (f"--filter-{input_type}", f"--filter_{input_type}")
        arglist_default = [
            f"{key}={value}"
            for (key, value) in pybids_inputs[input_type]["filters"].items()
        ]
        arglist_default_string = " ".join(arglist_default)

        filter_opts.add_argument(
            *argnames,
            nargs="+",
            action=KeyValue,
            help=f"(default: {arglist_default_string})",
        )

    # general parser for
    # --wildcards_{input_type} {wildcard1} {wildcard2} ...
    # create wildcards parsers, one for each input_type
    wildcards_opts = parser.add_argument_group(
        "INPUT WILDCARDS",
        "File path entities to use as wildcards in snakemake",
    )

    for input_type in pybids_inputs.keys():
        argnames = (f"--wildcards-{input_type}", f"--wildcards_{input_type}")
        arglist_default = [f"{wc}" for wc in pybids_inputs[input_type]["wildcards"]]
        arglist_default_string = " ".join(arglist_default)

        wildcards_opts.add_argument(
            *argnames,
            nargs="+",
            help=f"(default: {arglist_default_string})",
        )

    override_opts = parser.add_argument_group(
        "PATH OVERRIDE",
        (
            "Options for overriding BIDS by specifying absolute paths "
            "that include wildcards, e.g.: "
            "/path/to/my_data/{subject}/t1.nii.gz"
        ),
    )

    # create path override parser
    for input_type in pybids_inputs.keys():
        argnames = (f"--path-{input_type}", f"--path_{input_type}")
        override_opts.add_argument(*argnames, default=None)


def parse_snakebids_args(parser: argparse.ArgumentParser):
    all_args = parser.parse_known_args()
    return SnakebidsArgs(
        snakemake_args=all_args[1],
        # resolve all path items to get absolute paths
        args_dict={k: _resolve_path(v) for k, v in all_args[0].__dict__.items()},
        workflow_mode=all_args[0].workflow_mode,
        force=all_args[0].force_conversion,
        outputdir=Path(all_args[0].output_dir).resolve(),
        retrofit=all_args[0].retrofit,
    )


def _make_underscore_dash_aliases(name: str):
    """Generate --dash-arg aliases for --dash_args and vice versa

    If no dashes or underscores are in the argument name, a tuple containing just the
    original arg name will be returned

    Parameters
    ----------
    name : str
        Argument to generate conversion for

    Returns
    -------
    tuple of strings
        Converted args
    """
    match = re.match(r"^--(.+)$", name)
    if match:
        name_part = match.group(1)
        return {
            name,
            re.sub(r"\_", "-", name),
            f"--{re.sub(r'-', '_', name_part)}",
        }
    return {name}


def _resolve_path(path_candidate: Any) -> Any:
    """Helper function to resolve any paths or list
    of paths it's passed. Otherwise, returns the argument
    unchanged.

    Parameters
    ----------
    command : list, os.Pathlike, object
        command to run

    Returns
    -------
    list, Path, object
        If os.Pathlike or list  of os.Pathlike, the same paths resolved.
        Otherwise, the argument unchanged.
    """
    if isinstance(path_candidate, list):
        return [_resolve_path(p) for p in path_candidate]

    if isinstance(path_candidate, os.PathLike):
        return Path(path_candidate).resolve()

    return path_candidate
