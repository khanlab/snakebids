"""Tools to generate a Snakemake-based BIDS app.

This legacy module once had the core snakebids bidsapp implementation, but now is a
simple wrapper around :mod:`snakebids.bidsapp` with the
:class:`~snakebids.plugins.SnakemakeBidsApp` plugin. For new apps, this functionality
can be more flexibly implemented with:

.. code-block:: python

    from snakebids import bidsapp, plugins

    bidsapp.app([plugins.SnakemakeBidsApp(...)])

"""

from __future__ import annotations

import logging
import warnings
from os import PathLike
from pathlib import Path
from typing import Any, Callable

import attrs
import boutiques.creator as bc  # type: ignore

from snakebids import bidsapp
from snakebids import plugins as sb_plugins
from snakebids.bidsapp.run import _Runner

logger = logging.Logger(__name__)


@attrs.define
class SnakeBidsApp:
    """Snakebids app with config and arguments.

    Parameters
    ----------
    snakemake_dir : str | Path
        Root directory of the snakebids app, containing the config file and workflow
        files.
    plugins
        List of plugins to be registered.

        See :ref:`using-plugins` for more info.
    skip_parse_args
        DEPRECATED: no-op.
    parser
        DEPRECATED: no-op. (Historic: Parser including only the arguments
        specific to this Snakebids app, as specified in the config file. By
        default, it will use `create_parser()` from `cli.py`)
    configfile_path
        Relative path to config file (relative to snakemake_dir). By default,
        autocalculates based on snamake_dir
    snakefile_path
        Absolute path to the input Snakefile. By default, autocalculates based on
        snakemake_dir::

            join(snakemake_dir, snakefile_path)
    config
        DEPRECATED: no-op. (Historic: Contains all the configuration variables
        parsed from the config file and generated during the initialization of
        the SnakeBidsApp.)
    args
        DEPRECATED: no-op. (Historic: Arguments to use when running the app. By
        default, generated using the parser attribute, autopopulated with args
        from `config.py`)
    version
        DEPRECATED: no-op, use version plugin instead
    """

    snakemake_dir: Path = attrs.field(converter=Path)
    plugins: list[Callable[[SnakeBidsApp], None | SnakeBidsApp]] = attrs.Factory(list)
    skip_parse_args: bool = False
    _parser: Any = attrs.field(default=None, alias="parser")
    configfile_path: Path | None = attrs.field(
        default=None, converter=attrs.converters.optional(Path)
    )
    snakefile_path: Path | None = attrs.field(
        default=None, converter=attrs.converters.optional(Path)
    )
    _config: Any = attrs.field(default=None, alias="config")
    version: str | None = None
    args: Any = None

    _app_holder: _Runner | None = attrs.field(default=None, init=False)

    @property
    def _app(self):
        if self._app_holder is None:
            self._check_deprecations()

            self._app_holder = bidsapp.app(
                [sb_plugins.SnakemakeBidsApp(**self._get_args()), *self.plugins],
                description="Snakebids helps build BIDS Apps with Snakemake",
            )
        return self._app_holder

    def _check_deprecations(self):
        if self._parser is not None:
            msg = (
                "`SnakeBidsApp.parser` is deprecated and no longer has any effect. To "
                "modify the parser, use the new `bidsapp` module."
            )
            warnings.warn(msg, stacklevel=3)
        if self._config is not None:
            msg = (
                "`SnakeBidsApp.config` is deprecated and no longer has any effect. To "
                "modify the config, use the new `bidsapp` module."
            )
            warnings.warn(msg, stacklevel=3)
        if self.args is not None:
            msg = (
                "`SnakeBidsApp.args` is deprecated and no longer has any effect. To "
                "modify CLI arguments, use the new `bidsapp` module."
            )
            warnings.warn(msg, stacklevel=3)
        if self.version is not None:
            msg = (
                "`SnakeBidsApp.version` is deprecated and no longer has any effect. To "
                "explcitly set the app version, use the `snakebids.plugins.Version` "
                "plugin"
            )
            warnings.warn(msg, stacklevel=3)

    @property
    def config(self):
        """Get config dict (before arguments are parsed)."""
        return self._app.build_parser().config

    @property
    def parser(self):
        """Get parser."""
        return self._app.build_parser().parser

    def _get_args(self):
        args: dict[str, Any] = {}
        args["snakemake_dir"] = self.snakemake_dir
        if self.configfile_path is not None:
            args["configfile_path"] = self.configfile_path
        if self.snakefile_path is not None:
            args["snakefile_path"] = self.snakefile_path
        return args

    def run_snakemake(self) -> None:
        """Run snakemake with the given config, after applying plugins."""
        self._app.run()

    def create_descriptor(self, out_file: PathLike[str] | str) -> None:
        """Generate a boutiques descriptor for this Snakebids app."""
        new_descriptor = bc.CreateDescriptor(  # type: ignore
            self.parser, execname="run.py"
        )
        new_descriptor.save(out_file)  # type: ignore
