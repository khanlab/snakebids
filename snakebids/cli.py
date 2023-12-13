from __future__ import annotations

import argparse
import logging
import pathlib
import re
from collections.abc import Sequence
from typing import Any, Mapping, TypeVar, overload

import attr
import snakemake
from typing_extensions import override

from snakebids.exceptions import ConfigError, MisspecifiedCliFilterError
from snakebids.io.yaml import get_yaml_io
from snakebids.types import InputsConfig, OptionalFilter
from snakebids.utils.utils import to_resolved_path

# We define Path here in addition to pathlib to put both variables in globals()
# This way, users specifying a path type in their config.yaml can indicate
# either Path or pathlib.Path
Path = pathlib.Path

logger = logging.Logger(__name__)


class FilterParse(argparse.Action):
    """Class for parsing CLI filters in argparse."""

    # Constructor calling
    @override
    def __call__(
        self,
        parser: argparse.ArgumentParser,
        namespace: argparse.Namespace,
        values: str | Sequence[Any] | None,
        option_string: str | None = None,
    ):
        setattr(namespace, self.dest, {})
        if not values:
            return

        for pair in values:
            if "=" in pair:
                # split it into key and value
                key, value = pair.split("=", 1)
            elif ":" in pair:
                key, spec = pair.split(":", 1)
                spec = spec.lower()
                if spec == "optional":
                    value = OptionalFilter
                elif spec in ["required", "any"]:
                    value = True
                elif spec == "none":
                    value = False
                else:
                    # The flag isn't recognized
                    raise MisspecifiedCliFilterError(pair)
            else:
                raise MisspecifiedCliFilterError(pair)

            # assign into dictionary
            getattr(namespace, self.dest)[key] = value


class SnakemakeHelpAction(argparse.Action):
    """Class for printing snakemake usage in argparse."""

    @override
    def __call__(
        self,
        parser: argparse.ArgumentParser,
        namespace: argparse.Namespace,
        values: str | Sequence[Any] | None,
        option_string: str | None = None,
    ):
        snakemake.main(["-h"])  # type: ignore


@attr.frozen
class SnakebidsArgs:
    """Arguments for Snakebids App.

    Organizes the various options available for a generic Snakebids App, and store
    project specific arguments in a dict. Snakemake args are to be put in a list

    Attributes
    ----------
    force : bool
        Force output in a directory that already has contents
    outputdir : Path
        Directory to place outputs
    pybidsdb_dir : Path
        Directory to place pybids database
    snakemake_args : list of strings
        Arguments to pass on to Snakemake
    args_dict : dict[str, Any]
        Contains all the snakebids specific args. Meant to contain custom user args
        defined in config, as well as dynamic --filter-xx and --wildcard-xx args.
        These will eventually be printed in the new config.
    """

    force: bool
    outputdir: Path = attr.ib(converter=to_resolved_path)
    snakemake_args: list[str]
    args_dict: dict[str, Any]
    pybidsdb_dir: Path | None = None
    pybidsdb_reset: bool = False


def create_parser(include_snakemake: bool = False) -> argparse.ArgumentParser:
    """Generate basic Snakebids Parser.

    Includes the standard Snakebids arguments.
    """
    # The snakemake parser functionality did not seem to be implemented in the original
    # factoring. I left the logic here, but it should probably be removed if it's not
    # needed.
    if include_snakemake:
        # get snakemake parser
        smk_parser = snakemake.get_argument_parser()  # type: ignore

        # create parser
        parser = argparse.ArgumentParser(
            description="Snakebids helps build BIDS Apps with Snakemake",
            add_help=False,
            parents=[smk_parser],  # type: ignore
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
        help=argparse.SUPPRESS,
    )

    # We use -x as the alias because both -f and -F are taken by snakemake
    standard_group.add_argument(
        "--force-conversion",
        "--force_conversion",
        "-x",
        action="store_true",
        help=argparse.SUPPRESS,
    )

    standard_group.add_argument(
        "--pybidsdb-dir",
        "--pybidsdb_dir",
        action="store",
        help=(
            "Optional path to directory of SQLite databasefile for PyBIDS. "
            "If directory is passed and folder exists, indexing is skipped. "
            "If pybidsdb_reset is called, indexing will persist"
        ),
    )

    standard_group.add_argument(
        "--pybidsdb-reset",
        "--pybidsdb_reset",
        action="store_true",
        help=("Reindex existing PyBIDS SQLite database"),
    )

    # To be deprecated
    standard_group.add_argument(
        "--reset-db",
        "--reset_db",
        action="store_true",
        help=argparse.SUPPRESS,
    )

    standard_group.add_argument(
        "--force-output",
        "--force_output",
        action="store_true",
        help="Force output in a new directory that already has contents",
    )

    standard_group.add_argument(
        "--retrofit",
        action="store_true",
        help=argparse.SUPPRESS,
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


def _find_type(name: str, *, yamlsafe: bool = True) -> type[Any]:
    import importlib

    if name == "Path":
        return Path
    *module_name, obj_name = name.split(".") if "." in name else ("builtins", name)
    try:
        type_ = getattr(importlib.import_module(".".join(module_name)), obj_name)
    except (ImportError, AttributeError) as err:
        msg = f"{name} could not be resolved"
        raise ConfigError(msg) from err
    if not callable(type_):
        msg = f"{name} cannot be used as a type"
        raise ConfigError(msg)
    if yamlsafe and type_ not in get_yaml_io().representer.yaml_representers:
        msg = f"{name} cannot be serialized into yaml"
        raise ConfigError(msg)
    return type_


def add_dynamic_args(
    parser: argparse.ArgumentParser,
    parse_args: Mapping[str, Any],
    pybids_inputs: InputsConfig,
) -> None:
    """Add --filter-<comp> and --wildcards-<comp> arguments to CLI."""
    # create parser group for app options
    app_group = parser.add_argument_group("SNAKEBIDS", "Options for snakebids app")

    # update the parser with config options
    for name, arg in parse_args.items():
        # Convert type annotations from strings to class types
        # We first check that the type annotation is, in fact,
        # a str to allow the edge case where it's already
        # been converted
        if "type" in arg:
            arg_dict = {**arg, "type": _find_type(str(arg["type"]))}
        else:
            arg_dict = arg
        app_group.add_argument(
            *_make_underscore_dash_aliases(name),
            **arg_dict,  # type: ignore
        )

    # general parser for
    # --filter_{input_type} {key1}={value1} {key2}={value2}...
    # create filter parsers, one for each input_type
    filter_opts = parser.add_argument_group(
        "BIDS FILTERS",
        "Filters to customize PyBIDS get() as key=value pairs, or as "
        "key:{REQUIRED|OPTIONAL|NONE} (case-insensitive), to enforce the presence or "
        "absence of values for that key.",
    )

    for input_type in pybids_inputs:
        argnames = (f"--filter-{input_type}", f"--filter_{input_type}")
        filters = pybids_inputs[input_type].get("filters", {})
        arglist_default = [f"{key}={value}" for (key, value) in filters.items()]

        filter_opts.add_argument(
            *argnames,
            nargs="+",
            action=FilterParse,
            help=f"(default: {' '.join(arglist_default)})",
        )

    # general parser for
    # --wildcards_{input_type} {wildcard1} {wildcard2} ...
    # create wildcards parsers, one for each input_type
    wildcards_opts = parser.add_argument_group(
        "INPUT WILDCARDS",
        "File path entities to use as wildcards in snakemake",
    )

    for input_type in pybids_inputs:
        argnames = (f"--wildcards-{input_type}", f"--wildcards_{input_type}")
        arglist_default = [
            f"{wc}" for wc in pybids_inputs[input_type].get("wildcards", [])
        ]

        wildcards_opts.add_argument(
            *argnames,
            nargs="+",
            help=f"(default: {' '.join(arglist_default)})",
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
    for input_type in pybids_inputs:
        argnames = (f"--path-{input_type}", f"--path_{input_type}")
        override_opts.add_argument(*argnames, default=None)


def parse_snakebids_args(parser: argparse.ArgumentParser) -> SnakebidsArgs:
    """Parse built-in snakebids arguments."""
    all_args = parser.parse_known_args()
    if all_args[0].workflow_mode:
        logger.warning(
            "--workflow-mode is deprecated, and no longer has any effect. Workflow "
            "mode is automatially activated when choosing the snakemake root "
            "directory, or a subfolder of root/results, as the output directory."
        )
    if all_args[0].retrofit:
        logger.warning(
            "--retrofit is deprecated and no longer has any effect. To run this "
            "snakebids app on an old output (e.g. if snakebids raises an error because "
            "the directory already has contents, use the --force-output flag."
        )
    if all_args[0].force_conversion:
        logger.warning("--force-conversion is deprecated and no longer has any effect.")
    if all_args[0].reset_db:
        logger.warning(
            "--reset-db/--reset_db will be deprecated in a future release. To reset "
            "the pybids database, use the new --pybidsdb-reset flag instead."
        )
    return SnakebidsArgs(
        snakemake_args=all_args[1],
        # resolve all path items to get absolute paths
        args_dict={k: _resolve_path(v) for k, v in all_args[0].__dict__.items()},
        force=all_args[0].force_output,
        outputdir=Path(all_args[0].output_dir).resolve(),
        pybidsdb_dir=(
            None
            if all_args[0].pybidsdb_dir is None
            else Path(all_args[0].pybidsdb_dir).resolve()
        ),
        pybidsdb_reset=all_args[0].pybidsdb_reset or all_args[0].reset_db,
    )


def _make_underscore_dash_aliases(name: str) -> set[str]:
    """Generate --dash-arg aliases for --dash_args and vice versa.

    If no dashes or underscores are in the argument name, a tuple containing just the
    original arg name will be returned

    Parameters
    ----------
    name : str
        Argument to generate conversion for

    Returns
    -------
    set of strings
        Converted args
    """
    if match := re.match(r"^--(.+)$", name):
        name_part = match.group(1)
        return {
            name,
            re.sub(r"\_", "-", name),
            f"--{re.sub(r'-', '_', name_part)}",
        }
    return {name}


_T = TypeVar("_T")


@overload
def _resolve_path(path_candidate: Sequence[Any]) -> list[Any]:
    ...


@overload
def _resolve_path(path_candidate: _T) -> _T:
    ...


def _resolve_path(path_candidate: Any) -> Any:
    """Resolve paths or list of paths, or return argument unchanged.

    Parameters
    ----------
    command : list, os.Pathlike, object
        command to run

    Returns
    -------
        If os.Pathlike or list  of os.Pathlike, the same paths resolved.
        Otherwise, the argument unchanged.
    """
    if isinstance(path_candidate, Sequence) and not isinstance(path_candidate, str):
        return [
            _resolve_path(p)  # type: ignore[reportUnknownArgumentType]
            for p in path_candidate  # type: ignore[reportUnknownVariableType]
        ]

    if isinstance(path_candidate, Path):
        return Path(path_candidate).resolve()

    return path_candidate
