from __future__ import annotations

import argparse
import logging
from pathlib import Path
from typing import Any

import attrs

from snakebids import bidsapp
from snakebids.bidsapp.args import ArgumentGroups
from snakebids.plugins.base import PluginBase
from snakebids.utils.utils import DEPRECATION_FLAG

logger = logging.getLogger(__name__)


@attrs.define
class Pybidsdb(PluginBase):
    """Add CLI parameters to specify and reset a pybids database.

    Parameters
    ----------
    argument_group
        Specify title of the group to which arguments should be added


    CLI Arguments
    ~~~~~~~~~~~~~
    Two arguments are added to the CLI. These can be overriden by adding arguments
    with corresponding ``dests`` before this plugin is run:

    - ``plugins.pybidsdb.dir``: (:class:`~pathlib.Path`) Path of the database
    - ``plugins.pybidsdb.reset``: (:class:`bool`) Boolean indicating the database should
      be reset.

    After parsing, the above dests will be moved into ``config`` under the following
    names:

    - ``plugins.pybidsdb.dir`` → ``pybidsdb_dir``
    - ``plugins.pybidsdb.reset`` → ``pybidsdb_reset``

    This plugin only handles the CLI arguments, it does not do any actions with the
    database. The above config entries can be consumed by downstream processes.
    """

    argument_group: str | None = None

    PREFIX = "plugins.pybidsdb"

    def __eq__(self, other: Any):
        return isinstance(other, self.__class__)

    @bidsapp.hookimpl
    def add_cli_arguments(
        self, parser: argparse.ArgumentParser, argument_groups: ArgumentGroups
    ):
        """Add database parameters."""
        group = (
            argument_groups[self.argument_group]
            if self.argument_group is not None
            else parser
        )
        self.try_add_argument(
            group,
            "--pybidsdb-dir",
            "--pybidsdb_dir",
            action="store",
            type=Path,
            dest="dir",
            metavar="PATH",
            help=(
                "Optional path to directory of SQLite databasefile for PyBIDS. "
                "If directory is passed and folder exists, indexing is skipped. "
                "If pybidsdb_reset is called, indexing will persist"
            ),
        )

        self.try_add_argument(
            group,
            "--pybidsdb-reset",
            "--pybidsdb_reset",
            action="store_true",
            dest="reset",
            help="Reindex existing PyBIDS SQLite database",
        )

        # To be deprecated
        self.try_add_argument(
            group,
            "--reset-db",
            "--reset_db",
            action="store_true",
            dest="reset_db",
            help=argparse.SUPPRESS,
        )

    @bidsapp.hookimpl
    def update_cli_namespace(self, namespace: dict[str, Any], config: dict[str, Any]):
        """Assign database parameters to config."""
        # Update config with pybids settings
        reset = self.pop(namespace, "reset")
        reset_db = self.pop(namespace, "reset_db")
        if reset_db:
            logger.warning(
                "--reset-db/--reset_db will be deprecated in a future release. To "
                "reset the pybids database, use the new --pybidsdb-reset flag instead."
            )
        config["pybidsdb_reset"] = reset or reset_db
        pybidsdb_dir = self.pop(namespace, "dir")

        config["pybidsdb_dir"] = (
            None if pybidsdb_dir is None else pybidsdb_dir.resolve()
        )
        config["pybids_db_dir"] = (
            f"{DEPRECATION_FLAG}{config['pybidsdb_dir']}{DEPRECATION_FLAG}"
            if config["pybidsdb_dir"] is not None
            else None
        )
        config["pybids_db_reset"] = (
            f"{DEPRECATION_FLAG}{int(config['pybidsdb_reset'])}{DEPRECATION_FLAG}"
        )
