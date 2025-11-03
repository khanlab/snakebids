from __future__ import annotations

import itertools as it
from argparse import ArgumentParser
from collections.abc import Iterable
from pathlib import Path

import pytest
from hypothesis import given
from hypothesis import strategies as st

from snakebids.exceptions import ConfigError
from snakebids.plugins.bidsargs import BidsArgs
from tests import strategies as sb_st


class _St:
    positional_args = st.text().filter(lambda s: not s.startswith("-"))
    analysis_level_choices = positional_args | st.iterables(positional_args)


class TestAddCliArguments:
    @given(
        analysis_level_choices=_St.positional_args,
        analysis_level_choices_config=st.text(),
    )
    def test_analysis_level_attributes_are_mutually_exclusive(
        self,
        analysis_level_choices: str | Iterable[str],
        analysis_level_choices_config: str,
    ):
        with pytest.raises(ConfigError, match="cannot be simultaneously defined"):
            BidsArgs(
                analysis_level=True,
                analysis_level_choices=analysis_level_choices,
                analysis_level_choices_config=analysis_level_choices_config,
            ).add_cli_arguments(ArgumentParser(), {}, {})

    @given(analysis_level_choices=_St.positional_args)
    def test_analysis_levels_can_be_directly_defined(
        self, analysis_level_choices: str | Iterable[str]
    ):
        parser = ArgumentParser()
        a1, a2 = it.tee(analysis_level_choices)
        bidsargs = BidsArgs(analysis_level=True, analysis_level_choices=a1)
        bidsargs.add_cli_arguments(parser, {}, {})
        for arg in a2:
            nspc = parser.parse_args(["...", "...", arg])
            assert nspc.analysis_level == arg

    @given(
        analysis_level_choices=_St.positional_args,
        analysis_level_choices_config=st.text(),
    )
    def test_analysis_levels_can_be_defined_in_config(
        self,
        analysis_level_choices: str | Iterable[str],
        analysis_level_choices_config: str,
    ):
        parser = ArgumentParser()
        choices = list(analysis_level_choices)
        bidsargs = BidsArgs(
            analysis_level=True,
            analysis_level_choices_config=analysis_level_choices_config,
        )
        bidsargs.add_cli_arguments(parser, {}, {analysis_level_choices_config: choices})
        for arg in choices:
            nspc = parser.parse_args(["...", "...", arg])
            assert nspc.analysis_level == arg

    @given(
        derivatives=st.lists(
            st.text(sb_st.path_characters).filter(lambda s: not s.startswith("-")),
            min_size=1,
        )
    )
    def test_derivatives_converted_to_paths(self, derivatives: list[str]):
        parser = ArgumentParser()
        bidsargs = BidsArgs()
        bidsargs.add_cli_arguments(parser, {}, {})
        nspc = parser.parse_args(
            ["...", "...", "participant", "--derivatives", *derivatives]
        )
        assert nspc.derivatives == [Path(p) for p in derivatives]

    def test_derivatives_true_if_no_paths_given(self):
        parser = ArgumentParser()
        bidsargs = BidsArgs()
        bidsargs.add_cli_arguments(parser, {}, {})
        nspc = parser.parse_args(["...", "...", "participant", "--derivatives"])
        assert nspc.derivatives is True

    def test_derivatives_false_by_default(self):
        parser = ArgumentParser()
        bidsargs = BidsArgs()
        bidsargs.add_cli_arguments(parser, {}, {})
        nspc = parser.parse_args(["...", "...", "participant"])
        assert nspc.derivatives is False
