from __future__ import annotations

import argparse
from typing import Any

import pluggy

from snakebids.bidsapp.args import ArgumentGroups

hookspec = pluggy.HookspecMarker("snakebids")


@hookspec
def initialize_config(config: dict[str, Any]):
    """Modify config dict before creation of CLI parser.

    Config should be modified in place.

    Parameters
    ----------
    config
        Possibly empty dictionary of configuration values
    """


@hookspec
def add_cli_arguments(
    parser: argparse.ArgumentParser,
    config: dict[str, Any],
    argument_groups: ArgumentGroups,
):
    """Add any number of arguments to the argument parser.

    If special intervention is not made in the :func:`update_cli_namespace` hook,
    argument data will be automatically merged into ``config`` following argument
    parsing under the argument `dest`. To avoid naming conflicts, plugins should
    explicitly set the `dest` of each argument they add to
    ``plugins.<plugin_name>.<argument_name>``.

    Parameters
    ----------
    parser
        BIDS app CLI parser
    config
        Configuration dictionary preloaded with values from :func:`initialize_config`
    """


@hookspec(firstresult=True)
def get_argv(argv: list[str], config: dict[str, Any]) -> list[str] | None:
    """Set or modify the CLI parameters that will be parsed by the parser."""


@hookspec
def handle_unknown_args(args: list[str], config: dict[str, Any]):
    """If ``parse_known_args`` enabled, handle unknown arguments.

    Parameters
    ----------
    args
        List of unknown arguments parsed by the argparse parser
    config
        Configuration dictionary
    """


@hookspec
def update_cli_namespace(namespace: dict[str, Any], config: dict[str, Any]):
    """Interact with argument parsing results before they are merged into ``config``.

    The ``namespace`` contains the results from
    :meth:`~argparse.ArgumentParser.parse_args` (equivalent to
    :class:`vars(Namespace) <argparse.Namespace>`). Any modifications made to this
    :class:`dict` will be carried forward in app initialization. For instance, if
    an entry is deleleted from ``namespace``, it will not be available to downstream
    plugins or be copied into ``config``.

    Parameters
    ----------
    namespace
        the dictionary of values derived from the argparse namespace
    config
        Configuration dictionary
    """


@hookspec
def finalize_config(config: dict[str, Any]):
    """Perform modifications to the config following merging of parsed arguments.

    Parameters
    ----------
    config
        Configuration dictionary
    """


@hookspec
def run(config: dict[str, Any]):
    """Consume configuration and perform actions.

    This hook runs directly after :func:`finalize_config`. The config should not be
    modified.

    Parameters
    ----------
    config
        The finalized configuration dictionary
    """
