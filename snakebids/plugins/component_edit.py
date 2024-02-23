from __future__ import annotations

import argparse
import math
from pathlib import Path
from typing import Any, Sequence

import attrs
from typing_extensions import override

from snakebids import bidsapp
from snakebids.exceptions import MisspecifiedCliFilterError
from snakebids.plugins.base import PluginBase
from snakebids.types import OptionalFilter


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
            eq = pair.find("=")
            col = pair.find(":")
            delim = min(eq if eq >= 0 else math.inf, col if col >= 0 else math.inf)
            if delim is math.inf:
                raise MisspecifiedCliFilterError(pair)
            key = pair[:delim]
            value = pair[delim + 1 :]
            if delim == col:
                spec = value.lower()
                if spec == "optional":
                    value = OptionalFilter
                elif spec in ["required", "any"]:
                    value = True
                elif spec == "none":
                    value = False
                else:
                    # The flag isn't recognized
                    raise MisspecifiedCliFilterError(pair)

            # assign into dictionary
            getattr(namespace, self.dest)[key] = value


@attrs.define
class ComponentEdit(PluginBase):
    """Use CLI arguments to edit the filters, wildcards, and paths of components.

    Arguments are added based on the components specified in ``config``. For each
    component, a ``--filter-<comp_name>``, ``--wildcards-<comp_name>``, and a
    ``--path-<comp_name>`` argument will be added to the CLI. After parsing, these
    arguments are read and used to update the original component specification within
    config.

    CLI arguments created by this plugin cannot be overriden.

    Parameters
    ----------
    components_key
        Key of component specification within the config dictionary.
    """

    components_key: str = "pybids_inputs"

    PREFIX = "plugins.componentedit"

    def __eq__(self, other: Any):
        return isinstance(other, self.__class__)

    @bidsapp.hookimpl
    def add_cli_arguments(
        self, parser: argparse.ArgumentParser, config: dict[str, Any]
    ):
        """Add filter, wildcard, and path override arguments for each component."""
        if (pybids_inputs := config.get(self.components_key)) is None:
            return

        # general parser for
        # --filter_{input_type} {key1}={value1} {key2}={value2}...
        # create filter parsers, one for each input_type
        filter_opts = parser.add_argument_group(
            "BIDS FILTERS",
            "Filters to customize PyBIDS get() as key=value pairs, or as "
            "key:{REQUIRED|OPTIONAL|NONE} (case-insensitive), to enforce the presence "
            "or absence of values for that key.",
        )

        for input_type in pybids_inputs:
            argnames = (f"--filter-{input_type}", f"--filter_{input_type}")
            filters = pybids_inputs[input_type].get("filters", {})
            arglist_default = [f"{key}={value}" for (key, value) in filters.items()]

            self.add_argument(
                filter_opts,
                *argnames,
                nargs="+",
                action=FilterParse,
                dest=f"{self.PREFIX}.filter.{input_type}",
                metavar="ENTITY=VALUE",
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

            self.add_argument(
                wildcards_opts,
                *argnames,
                nargs="+",
                dest=f"{self.PREFIX}.wildcards.{input_type}",
                metavar="WILDCARD",
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
            self.add_argument(
                override_opts,
                *argnames,
                default=None,
                metavar="PATH",
                dest=f"{self.PREFIX}.path.{input_type}",
            )

    @bidsapp.hookimpl
    def update_cli_namespace(self, namespace: dict[str, Any], config: dict[str, Any]):
        """Apply provided overrides to the component configuration."""
        # argparse adds filter_{input_type}
        # we want to update the pybids_inputs dict with this, then remove the
        # filter_{input_type} dict
        if (pybids_inputs := config.get(self.components_key)) is None:
            return
        for input_type in pybids_inputs:
            arg_filter_dict = self.pop(namespace, f"filter.{input_type}", None)
            if arg_filter_dict is not None:
                pybids_inputs[input_type].setdefault("filters", {})
                for entity, filter_ in arg_filter_dict.items():
                    if filter_ is OptionalFilter:
                        pybids_inputs[input_type]["filters"].pop(entity, None)
                    else:
                        pybids_inputs[input_type]["filters"][entity] = filter_

        # add cmdline defined wildcards from the list:
        # wildcards_{input_type}
        for input_type in pybids_inputs:
            wildcards_list = self.pop(namespace, f"wildcards.{input_type}", None)
            if wildcards_list is not None:
                pybids_inputs[input_type].setdefault("wildcards", [])
                pybids_inputs[input_type]["wildcards"] += wildcards_list

        # add custom input paths to
        # config['pybids_inputs'][input_type]['custom_path']
        for input_type in pybids_inputs:
            custom_path = self.pop(namespace, f"path.{input_type}", None)
            if custom_path is not None:
                pybids_inputs[input_type]["custom_path"] = Path(custom_path).resolve()
