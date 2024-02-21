from __future__ import annotations

import argparse
import re
from pathlib import Path
from typing import Any

import attrs

from snakebids import bidsapp
from snakebids.exceptions import ConfigError
from snakebids.io.yaml import get_yaml_io


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


@attrs.define
class CliConfig:
    """Configure CLI arguments directly in config.

    Arguments are provided in config in a dictionary stored under ``cli_key``. Each
    entry maps the name of the argument to a set of valid arguments for
    :meth:`~argparse.ArgumentParser.add_argument()`.

    This plugin will attempt to be the first to add arguments, and thus can be used
    to override arguments from other compatible plugins, such as
    :class:`~snakebids.plugins.bidsargs.BidsArgs`

    Parameters
    ----------
    cli_key
        Key of dict containing arguments

    Example
    -------

    .. code-block:: yaml

        parse_args:
            --tunable-parameter:
                help: |
                    A parameter important to the analysis that you can be set from the
                    commandline. If not set, a sensible default will be used
                default: 5
                type: float
            --alternate-mode:
                help: |
                    A flag activating a secondary feature of the workflow
                action: store_true
            --derivatives:
                help: |
                    An alternate help message for --derivatives, which will override
                    the help message from argument given by ``BidsArgs``. Note that we
                    must again specify the ``nargs`` and ``type`` for the param
                type: Path
                nargs: "*"
    """

    cli_key: str = "parse_args"

    def __eq__(self, other: Any):
        return isinstance(other, self.__class__)

    @bidsapp.hookimpl(tryfirst=True)
    def add_cli_arguments(
        self, parser: argparse.ArgumentParser, config: dict[str, Any]
    ):
        """Add arguments from config."""
        if (parse_args := config.get(self.cli_key)) is None:
            return

        # update the parser with config options
        for name, arg in parse_args.items():
            import pathlib as pathlib  # noqa: PLC0414
            from pathlib import Path as Path  # noqa: PLC0414

            # Convert type annotations from strings to class types
            # We first check that the type annotation is, in fact,
            # a str to allow the edge case where it's already
            # been converted
            if "type" in arg:
                arg_dict = {**arg, "type": _find_type(str(arg["type"]))}
            else:
                arg_dict = arg
            parser.add_argument(
                *_make_underscore_dash_aliases(name),
                **arg_dict,  # type: ignore
            )
