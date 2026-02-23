from __future__ import annotations

import argparse
from collections.abc import Sequence
from pathlib import Path
from typing import Any, Final, cast

import attrs
from typing_extensions import override

from snakebids import bidsapp
from snakebids.plugins.base import PluginBase
from snakebids.types import OptionalFilter
from snakebids.utils.utils import text_fold


class FilterParseError(Exception):
    """Raised when a magic CLI filter cannot be parsed."""

    def __init__(self, filter: str, msg: str):
        super().__init__(
            "The following filter provided by the CLI could not be parsed: "
            f"'{filter}'.\n\t{msg}".expandtabs(4)
        )

    @classmethod
    def invalid_spec(cls, filter: str, spec: str):
        """Raise if spec not recognized."""
        return cls(
            filter,
            f"':{spec}' is not a valid filter method. Must be one of 'required', "
            "'optional', 'none', 'match', or 'search'.",
        )

    @classmethod
    def missing_value(cls, filter: str, key: str, spec: str):
        """Raise if no value provided."""
        return cls(
            filter, f"':{spec}' requires a value, specified as '{key}:{spec}=VALUE'."
        )

    @classmethod
    def only_key(cls, filter: str, key: str):
        """Raise if only a key provided."""
        return cls(filter, "Filters must be specified as ENTITY[:METHOD]=VALUE.")

    @classmethod
    def unneeded_value(cls, filter: str, key: str, spec: str, value: str):
        """Raise if value provided to a boolean filter method."""
        return cls(
            filter, f"'{key}:{spec}' should not be given a value (got '={value}')"
        )


class FilterParse(argparse.Action):
    """Class for parsing CLI filters in argparse."""

    boolean_filters: Final = {
        "optional": OptionalFilter,
        "required": True,
        "any": True,
        "none": False,
    }

    # Constructor calling
    @override
    def __call__(
        self,
        parser: argparse.ArgumentParser,
        namespace: argparse.Namespace,
        values: str | Sequence[Any] | None,
        option_string: str | None = None,
    ):
        if not values:
            return

        result: dict[str, Any] = getattr(namespace, self.dest, None) or {}
        for pair in values:
            if "=" in pair:
                # split it into key and value
                key, value = cast("tuple[str, str]", pair.split("=", 1))
            else:
                key = pair
                value = None
            if ":" in key:
                key, spec = cast("tuple[str, str]", key.split(":", 1))
                spec = spec.lower()
                if spec in self.boolean_filters:
                    if value is not None:
                        raise FilterParseError.unneeded_value(pair, key, spec, value)
                    value = self.boolean_filters[spec]
                elif spec in {"match", "search"}:
                    if value is None:
                        raise FilterParseError.missing_value(pair, key, spec)
                    value = {spec: value}
                else:
                    raise FilterParseError.invalid_spec(pair, spec)
            if value is None:
                raise FilterParseError.only_key(pair, key)

            # assign into dictionary
            result[key] = value

        setattr(namespace, self.dest, result)


@attrs.define
class ComponentEdit(PluginBase):
    """Use CLI arguments to edit the filters, wildcards, and paths of components.

    Arguments are added based on the components specified in ``config``. For each
    component, a ``--filter-<comp_name>``, ``--wildcards-<comp_name>``, and a
    ``--path-<comp_name>`` argument will be added to the CLI. After parsing, these
    arguments are read and used to update the original component specification within
    config.

    Filters are specified on the CLI using ``ENTITY[:METHOD][=VALUE]``, as follows:

    1. ``ENTITY=VALUE`` selects paths based on an exact value match.
    2. ``ENTITY:match=REGEX`` and ``ENTITY:search=REGEX`` selects paths using regex
       with :func:`re.match` and :func:`re.search` respectively. This syntax can be used
       to select multiple values (e.g. ``'session:match=01|02'``).
    3. ``ENTITY:required`` selects all paths with the entity, regardless of value.
    4. ``ENTITY:none`` selects all paths without the entity.
    5. ``ENTITY:any`` removes filters for the entity.

    CLI arguments created by this plugin cannot be overridden.

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
            text_fold(
                """
                Update filters for input components. Each filter can be specified as a
                ENTITY=VALUE pair to select an value directly. To use regex filtering,
                ENTITY:match=REGEX or ENTITY:search=REGEX can be used for re.match() or
                re.search() respectively. Regex can also be used to select multiple
                values, e.g. 'session:match=01|02'. ENTITY:required and ENTITY:none can
                be used to require or prohibit presence of an entity in selected paths,
                respectively. ENTITY:optional can be used to remove a filter.
                """
            ),
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
                metavar="ENTITY[:METHOD][=VALUE]",
                help=f"(default: {' '.join(arglist_default)})",
            )

        # general parser for
        # --wildcards_{input_type} {wildcard1} {wildcard2} ...
        # create wildcards parsers, one for each input_type
        wildcards_opts = parser.add_argument_group(
            "INPUT WILDCARDS",
            "Provide entities to be used as wildcards.",
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
