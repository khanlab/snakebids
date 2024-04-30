from __future__ import annotations

import argparse
import importlib.metadata as impm
import logging
import sys
from pathlib import Path
from typing import Any, Callable, Sequence, TypeVar

import attrs
import more_itertools as itx
from typing_extensions import overload, override

from snakebids import bidsapp
from snakebids.exceptions import ConfigError, RunError
from snakebids.io.config import write_config
from snakebids.plugins.bidsargs import BidsArgs
from snakebids.plugins.cli_config import CliConfig
from snakebids.plugins.component_edit import ComponentEdit
from snakebids.plugins.pybidsdb import Pybidsdb
from snakebids.snakemake_compat import load_configfile
from snakebids.snakemake_compat import main as snakemake_main
from snakebids.utils.output import (
    prepare_bidsapp_output,
    write_output_mode,
)
from snakebids.utils.utils import to_resolved_path

__all__ = ["SnakemakeBidsApp"]

logger = logging.getLogger(__name__)


SNAKEFILE_CHOICES = [
    "Snakefile",
    "snakefile",
    "workflow/Snakefile",
    "workflow/snakefile",
]


CONFIGFILE_CHOICES = [
    "config/snakebids.yml",
    "config/snakebids.yaml",
    "config/snakebids.json",
    "snakebids.yml",
    "snakebids.yaml",
    "snakebids.json",
    "config.yml",
    "config.yaml",
    "config.json",
    "config/config.json",
    "config/config.yml",
    "config/config.yaml",
]


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
        snakemake_main(["-h"])  # type: ignore


def _get_file_paths(
    choices: list[str], file_name: str
) -> Callable[[SnakemakeBidsApp], Path]:
    def wrapper(self: SnakemakeBidsApp):
        for path in choices:
            if (self.snakemake_dir / path).exists():
                return Path(self.snakemake_dir, path)

        msg = f"Error: no {file_name} file found, tried {', '.join(choices)}."
        raise ConfigError(msg)

    return wrapper


_T = TypeVar("_T")


@overload
def _resolve_path(path_candidate: Sequence[Any]) -> list[Any]: ...


@overload
def _resolve_path(path_candidate: _T) -> _T: ...


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


@attrs.define
class SnakemakeBidsApp:
    """Snakebids app with config and arguments.

    .. currentmodule:: snakebids.plugins

    Loads the :class:`~bidsargs.BidsArgs`, :class:`~cli_config.CliConfig`,
    :class:`~pybidsdb.Pybidsdb`, and :class:`~component_edit.ComponentEdit` plugins
    as dependencies.

    Parameters
    ----------
    snakemake_dir : str | Pathlike[str]
        Root directory of the snakebids app, containing the config file and workflow
        files.
    configfile_path
        Path of config file. By default, looks for ``config.yaml``, ``config.json``,
        ``snakebids.yaml``, or ``snakebids.json`` either within the ``snakemake_dir`` or
        within ``snakemake_dir/config``
    configfile_outpath
        Path of output configfile that will be consumed by the snakemake workflow. This
        is only necessary if the input configfile is not within the snakemake
        directory. If given, it should be the same as the path specified under
        ``configfile: ...` within the ``Snakefile``
    snakefile_path
        Absolute path to the input Snakefile. By default, looks for ``Snakefile`` within
        the ``snakemake_dir`` or ``snakemake_dir/workflow``
    """

    snakemake_dir: Path = attrs.field(converter=to_resolved_path)
    configfile_path: Path | None = attrs.field(
        default=attrs.Factory(
            _get_file_paths(CONFIGFILE_CHOICES, "config"), takes_self=True
        ),
        kw_only=True,
    )
    configfile_outpath: Path | None = attrs.field(default=None, kw_only=True)
    snakefile_path: Path = attrs.field(
        default=attrs.Factory(
            _get_file_paths(SNAKEFILE_CHOICES, "Snakefile"), takes_self=True
        ),
        kw_only=True,
    )

    # use lambda so that it can be patched by fakefs in tests
    cwd: Path = attrs.field(factory=lambda: Path(), init=False)
    """The current working directory for running snakemake.

    This is where the ``.snakemake`` metadata folder will be placed.
    """

    force_output: bool = attrs.field(default=False, init=False)
    """Allow specifying outputs in unrecognized, non-empty directories.

    After the first run, a ``.snakebids`` file will be created in the directory to mark
    the directory as safe for future runs.
    """

    DEPENDENCIES = (
        Pybidsdb(),
        ComponentEdit(),
        CliConfig(),
        BidsArgs(),
    )

    @classmethod
    def create_empty(cls):
        """Create empty instance of plugin."""
        return cls(
            snakemake_dir=Path(),
            configfile_path=None,
            snakefile_path=Path("Snakefile"),
        )

    @bidsapp.hookimpl
    def initialize_config(self, config: dict[str, Any]):
        """Read config file and load into config dict."""
        if self.configfile_path is not None:
            config.update(
                load_configfile(str(self.snakemake_dir / self.configfile_path))
            )

    @bidsapp.hookimpl
    def add_cli_arguments(self, parser: argparse.ArgumentParser):
        """Add snakemake help and force_output arguments."""
        parser.add_argument(
            "--help-snakemake",
            "--help_snakemake",
            nargs=0,
            action=SnakemakeHelpAction,
            help=(
                "Options to Snakemake can also be passed directly at the "
                "command-line, use this to print Snakemake usage"
            ),
        )
        parser.add_argument(
            "--force-output",
            "--force_output",
            action="store_true",
            help="Force output in a new directory that already has contents",
        )

    @bidsapp.hookimpl
    def handle_unknown_args(self, args: list[str], config: dict[str, Any]):
        """Add snakemake args to config."""
        config["snakemake_args"] = args

    @bidsapp.hookimpl
    def update_cli_namespace(self, namespace: dict[str, Any], config: dict[str, Any]):
        """Resolve Paths within the namespace and set target."""
        if "output_dir" in namespace:
            namespace["output_dir"] = Path(namespace["output_dir"])
        if "bids_dir" in namespace:
            namespace["bids_dir"] = Path(namespace["bids_dir"])
        for key in namespace:
            namespace[key] = _resolve_path(namespace[key])
        namespace.pop("help_snakemake", None)
        self.force_output = namespace.pop("force_output", False)
        if "targets_by_analysis_level" in config and "analysis_level" in namespace:
            config["snakemake_target"] = config["targets_by_analysis_level"][
                namespace["analysis_level"]
            ]
        else:
            config["snakemake_target"] = None

    @bidsapp.hookimpl
    def finalize_config(self, config: dict[str, Any]):
        """Perform final steps for snakemake workflows.

        Expects to find ``output_dir`` in config as a fully resolved
        :class:`~pathlib.Path`

        - Modify ``output_dir`` to `output/results` if output is set to the app dir
        - Set ``self.cwd``
        - Write config file into output dir
        - Add snakemake specific parameters to config

        """
        # First, handle outputs in snakebids_root or results folder
        try:
            # py3.9 has the Path.is_relative() function. But as long as we support py38
            # and lower, this is the easiest way
            relative = config["output_dir"].relative_to(self.snakemake_dir / "results")
        except ValueError:
            relative = None

        if self.snakemake_dir == config["output_dir"] or relative is not None:
            write_output_mode(self.snakemake_dir / ".snakebids", "workflow")

            self.cwd = self.snakemake_dir
            root = Path("results", relative or "")

            if config["output_dir"] == self.snakemake_dir.resolve():
                config["output_dir"] /= "results"
                # Print a friendly warning that the output directory will change
                logger.info(
                    "You specified your output to be in the snakebids directory, so "
                    "we're automatically putting your outputs in the results "
                    "subdirectory.\nYou'll find your results in `%s`",
                    (self.snakemake_dir / "results").resolve(),
                )

        # else, we run in bidsapp mode
        else:
            try:
                prepare_bidsapp_output(config["output_dir"], self.force_output)
            except RunError as err:
                print(err.msg, file=sys.stderr)
                sys.exit(1)
            self.cwd = config["output_dir"]
            root = Path()

        configfile_path = self.configfile_path or self.snakemake_dir / "snakebids.yaml"
        if self.configfile_outpath is None:
            try:
                self.configfile_outpath = (
                    self.cwd
                    / configfile_path.resolve().relative_to(
                        self.snakemake_dir.resolve()
                    )
                )
            except ValueError:
                self.configfile_outpath = configfile_path

        version = config.get("plugins.version.version")
        if version == "unknown":
            logger.warning(
                "App version could not be identified. This usually occurs because the "
                "app was not installed with pip. It will be recorded in output as "
                "'unknown'."
            )
        # Write the config file
        write_config(
            config_file=self.configfile_outpath,
            data=dict(
                config,
                snakemake_version=impm.version("snakemake"),
                snakebids_version=impm.version("snakebids"),
                root=root,
                app_version=version or "unknown",
                snakemake_dir=self.snakemake_dir,
                snakefile=self.snakefile_path,
            ),
            force_overwrite=True,
        )

    @bidsapp.hookimpl
    def run(self, config: dict[str, Any]):
        """Run snakemake with the given config, after applying plugins."""
        snakemake_main(  # type: ignore
            [
                *filter(
                    None,
                    [
                        *itx.always_iterable(config["snakemake_target"]),
                        "--snakefile",
                        str(self.snakefile_path),
                        "--directory",
                        str(self.cwd),
                        "--configfile",
                        str(self.configfile_outpath),
                        *config["snakemake_args"],
                    ],
                )
            ]
        )
