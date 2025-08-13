from __future__ import annotations

import argparse
from pathlib import Path
from typing import Literal

import more_itertools as itx
import pytest
from hypothesis import given
from hypothesis import strategies as st

from snakebids.exceptions import ConfigError
from snakebids.plugins.cli_config import CliConfig


class TestMakesUnderscoreDashAliases:
    @given(
        args=st.lists(
            st.from_regex(r"[^-_\n]+-[^-_\n]+_[^-_\n]+", fullmatch=True),
            max_size=5,
            unique=True,
        )
    )
    @pytest.mark.parametrize("hot", [("dashed",), ("scored",)])
    def test_underscores_to_dashes(
        self, args: list[str], hot: Literal["dashed", "scored"]
    ):
        cli = CliConfig()
        parser = argparse.ArgumentParser()
        dashed = ["--" + s.replace("_", "-") for s in args]
        scored = ["--" + s.replace("-", "_") for s in args]
        curr = dashed if hot == "dashed" else scored
        opp = scored if hot == "dashed" else dashed
        config = {"parse_args": {arg: {"action": "store_true"} for arg in curr}}
        cli.add_cli_arguments(parser, config)
        nm = parser.parse_args(opp)
        assert all(nm.__dict__.values())
        assert len(nm.__dict__) == len(args)

    @given(arg=st.text().filter(lambda s: not s.startswith("-") and len(s)))
    def test_positional_args_added_without_conversion(self, arg: str):
        cli = CliConfig()
        parser = argparse.ArgumentParser()
        config = {"parse_args": {arg: {"help": "..."}}}
        cli.add_cli_arguments(parser, config)
        nm = parser.parse_args(["..."])
        assert itx.one(nm.__dict__.values()) == "..."

    def test_fails_if_undefined_type_given(self):
        cli = CliConfig()
        config = {"parse_args": {"--new-param": {"type": "UnheardClass"}}}
        parser = argparse.ArgumentParser()
        with pytest.raises(ConfigError, match="could not be resolved"):
            cli.add_cli_arguments(parser, config)


def test_works_when_key_not_present():
    cli = CliConfig()
    parser = argparse.ArgumentParser()
    cli.add_cli_arguments(parser, {})
    assert len(parser._actions) == 1


def test_convert_arg_to_builtin():
    cli = CliConfig()
    parser = argparse.ArgumentParser()
    new_args = {
        "--new-param": {
            "help": "Generic Help Message",
            "type": "int",
        }
    }
    cli.add_cli_arguments(parser, {"parse_args": new_args})
    args = parser.parse_args(["--new-param", "12"])
    assert isinstance(args.new_param, int)


def test_path_works_as_builtin():
    cli = CliConfig()
    parser = argparse.ArgumentParser()
    new_args = {
        "--new-param": {
            "help": "Generic Help Message",
            "type": "Path",
        }
    }
    cli.add_cli_arguments(parser, {"parse_args": new_args})
    args = parser.parse_args(["--new-param", "foo"])
    assert args.new_param == Path("foo")


def test_non_serialiable_type_raises_error():
    cli = CliConfig()
    parser = argparse.ArgumentParser()
    new_args = {
        "--new-param": {
            "help": "Generic Help Message",
            "type": "snakebids.utils.utils.BidsEntity",
        }
    }
    with pytest.raises(ConfigError, match="cannot be serialized into yaml"):
        cli.add_cli_arguments(parser, {"parse_args": new_args})


def test_using_module_as_type_gives_error():
    cli = CliConfig()
    parser = argparse.ArgumentParser()
    new_args = {
        "--new-param": {
            "help": "Generic Help Message",
            "type": "snakebids.utils.utils",
        }
    }
    with pytest.raises(ConfigError, match="cannot be used as a type"):
        cli.add_cli_arguments(parser, {"parse_args": new_args})


def test_using_class_method_as_type_gives_error():
    cli = CliConfig()
    parser = argparse.ArgumentParser()
    new_args = {
        "--new-param": {
            "help": "Generic Help Message",
            "type": "snakebids.utils.utils.BidsEntity.from_tag",
        }
    }
    with pytest.raises(ConfigError, match="could not be resolved"):
        cli.add_cli_arguments(parser, {"parse_args": new_args})
