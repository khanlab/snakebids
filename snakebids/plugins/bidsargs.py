from __future__ import annotations

import argparse
from pathlib import Path
from typing import Any, Iterable, Sequence

import attrs
import more_itertools as itx

from snakebids import bidsapp
from snakebids.bidsapp.args import ArgumentGroups
from snakebids.exceptions import ConfigError
from snakebids.plugins.base import PluginBase


def _list_or_none(arg: str | Iterable[str] | None) -> list[str] | None:
    return arg if arg is None else list(itx.always_iterable(arg))


class _Derivative(argparse.Action):
    def __call__(
        self,
        parser: argparse.ArgumentParser,
        namespace: argparse.Namespace,
        values: str | Sequence[Any] | None,
        option_string: str | None = None,
    ):
        # the below two conditions are to cover the typing. They are not actually
        # encountered at this time in the code
        if values is None:  # pragma: no cover
            result = False
        elif isinstance(values, str):  # pragma: no cover
            result = [Path(values)]
        elif not len(values):
            result = True
        else:
            result = [Path(p) for p in values]
        setattr(namespace, self.dest, result)


@attrs.define
class BidsArgs(PluginBase):
    """Add basic BIDSApp arguments.

    Parameters
    ----------
    argument_group
        Specify title of the group to which arguments should be added
    bids_dir
        Indicate if bids_dir (first bids argument) should be defined
    output_dir
        Indicate if output_dir (second bids argument) should be defined
    analysis_level
        Indicate if output_dir (third bids argument) should be defined
    analysis_level_choices
        List of valid analysis levels
    analysis_level_choices_config
        Configuration key containing analysis level choices
    participant_label
        Indicate if ``participant_label`` should be defined. Used to filter specific
        subjects for processesing
    exclude_participant_label
        Indicate if ``exclude_participant_label`` should be defined. Used to excluded
        specific subjects from processesing
    derivatives
        Indicate if ``derivatives`` should be defined. Used to allow automatic
        derivative indexing or specify paths to derivatives.


    CLI Arguments
    ~~~~~~~~~~~~~
    All arguments are added by default, but can be disabled using their respective
    parameters. Additinally, arguments can be directly overriden by adding arguments
    to the parser before the plugin runs, using the following ``dests``:

    - ``bids_dir``: The input bids directory
    - ``output_dir``: The output bids directory
    - ``analysis_level``: The level of analysis to perform (usually ``participant`` or
      ``group``)
    - ``participant_label``: Collection of subject labels to include in analysis
    - ``exclude_participant_label``: Collection of subject labels to exclude from
      analysis
    - ``derivatives``: Collection of derivative folder paths to include in bids
      indexing.

    Note that only the ``dest`` field matters when overriding arguments, all other
    feature, including positional versus optional, are incidental. Other than adding
    arguments for the above ``dests``, ``BidsArgs`` makes no assumptions about these
    arguments. Other plugins, however, may expect to find data in the parser namespace
    or config corresponding to the ``dest`` names.

    .. warning::

        Overriding just one or two of the positional arguments may alter the order,
        preventing the app from being called correctly. Thus, if any of the positional
        args are being overriden, they all should be.


    Analysis levels
    ~~~~~~~~~~~~~~~
    Valid analysis levels can be set using the ``analysis_level_choices`` key or the
    ``analysis_level_choices_config`` key. If both are set, an error will be raised. If
    neither are set, the config will first be searched for the ``analysis_levels`` key.
    If this is not found, analysis levels will be set to ``["participant", "group"]``.
    """

    argument_group: str = ""
    bids_dir: bool = True
    output_dir: bool = True
    analysis_level: bool = True
    analysis_level_choices: list[str] | None = attrs.field(
        default=None, converter=_list_or_none
    )
    analysis_level_choices_config: str | None = None
    participant_label: bool = True
    exclude_participant_label: bool = True
    derivatives: bool = True

    def __eq__(self, other: Any):
        return isinstance(other, self.__class__)

    @bidsapp.hookimpl
    def add_cli_arguments(
        self,
        parser: argparse.ArgumentParser,
        argument_groups: ArgumentGroups,
        config: dict[str, Any],
    ):
        """Add arguments from config."""
        group = argument_groups[self.argument_group] if self.argument_group else parser

        if self.bids_dir:
            self.try_add_argument(
                group,
                dest="bids_dir",
                help="The directory with the input dataset formatted according "
                "to the BIDS standard.",
                type=Path,
            )
        if self.output_dir:
            self.try_add_argument(
                group,
                dest="output_dir",
                help="The directory where the output files "
                "should be stored. If you are running group level analysis "
                "this folder should be prepopulated with the results of the "
                "participant level analysis.",
                type=Path,
            )
        if self.analysis_level:
            if (
                self.analysis_level_choices is not None
                and self.analysis_level_choices_config is not None
            ):
                msg = (
                    "`analysis_level_choices` and `analysis_level_choices_config` "
                    "cannot be simultaneously defined"
                )
                raise ConfigError(msg)
            if self.analysis_level_choices is not None:
                analysis_levels = self.analysis_level_choices
            elif self.analysis_level_choices_config is not None:
                analysis_levels = config.get(self.analysis_level_choices_config)
            elif (analysis_levels := config.get("analysis_levels")) is None:
                analysis_levels = ["participant", "group"]
            self.try_add_argument(
                group,
                dest="analysis_level",
                help="Level of the analysis that will be performed.",
                choices=analysis_levels,
            )
        if self.participant_label:
            self.try_add_argument(
                group,
                "--participant-label",
                "--participant_label",
                help="The label(s) of the participant(s) that should be analyzed. The "
                "label corresponds to sub-<participant_label> from the BIDS sec "
                '(so it does not include "sub-"). If this parameter is not '
                "provided all subjects should be analyzed. Multiple "
                "participants can be specified with a space separated list.",
                dest="participant_label",
                metavar="LABEL",
                nargs="+",
            )
        if self.exclude_participant_label:
            self.try_add_argument(
                group,
                "--exclude-participant-label",
                "--exclude_participant_label",
                help="The label(s) of the participant(s) that should be excluded. The "
                "label corresponds to sub-<participant_label> from the BIDS spec "
                '(so it does not include "sub-"). If this parameter is not '
                "provided all subjects should be analyzed. Multiple "
                "participants can be specified with a space separated list.",
                metavar="LABEL",
                dest="exclude_participant_label",
                nargs="+",
            )
        if self.derivatives:
            self.try_add_argument(
                group,
                "--derivatives",
                help="Path(s) to a derivatives dataset, for folder(s) that contains "
                "multiple derivatives datasets",
                nargs="*",
                dest="derivatives",
                metavar="PATH",
                action=_Derivative,
                default=False,
            )
