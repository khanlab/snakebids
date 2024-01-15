from __future__ import annotations

import copy
import itertools as it
import re
import sys
from argparse import ArgumentParser, Namespace
from collections.abc import Sequence
from pathlib import Path
from typing import ClassVar, Mapping

import hypothesis.strategies as st
import pytest
from hypothesis import given
from pytest_mock.plugin import MockerFixture

from snakebids.cli import (
    _resolve_path,
    add_dynamic_args,
    create_parser,
    parse_snakebids_args,
)
from snakebids.exceptions import ConfigError
from snakebids.tests import strategies as sb_st
from snakebids.tests.helpers import allow_function_scoped
from snakebids.types import InputsConfig, OptionalFilter

from .mock.config import parse_args, pybids_inputs


@pytest.fixture
def parser():
    p = create_parser()
    add_dynamic_args(p, copy.deepcopy(parse_args), copy.deepcopy(pybids_inputs))
    return p


class TestResolvePath:
    @pytest.fixture
    def arg_dict(self) -> dict[str, str | list[str]]:
        return {
            "bids_dir": "path/to/input",
            "output_dir": "path/to/output",
            "analysis_level": "participant",
            "--derivatives": ["path/to/deriv1", "path/to/deriv2"],
        }

    def test_does_not_change_dict_without_paths(
        self, arg_dict: Mapping[str, str | Sequence[str]]
    ):
        arg_dict_copy = copy.deepcopy(arg_dict)
        resolved = {key: _resolve_path(value) for key, value in arg_dict.items()}
        assert resolved == arg_dict_copy

    def test_resolves_all_paths(
        self, arg_dict: dict[str, str | Path | Sequence[str | Path]]
    ):
        derivative_paths = [Path("path/to/deriv1"), Path("path/to/deriv2")]
        arg_dict["--derivatives"] = derivative_paths
        arg_dict_copy = copy.deepcopy(arg_dict)
        arg_dict_copy["--derivatives"] = [p.resolve() for p in derivative_paths]
        resolved = {key: _resolve_path(value) for key, value in arg_dict.items()}
        assert resolved == arg_dict_copy

    def test_filter_dict(self, arg_dict: dict[str, dict[str, str]]):
        filter_dict = {"test_key": "test_value"}
        arg_dict["filter_test"] = filter_dict
        arg_dict_copy = copy.deepcopy(arg_dict)
        resolved = {key: _resolve_path(value) for key, value in arg_dict.items()}
        assert resolved == arg_dict_copy


class TestAddDynamicArgs:
    mock_args_special: ClassVar[list[str]] = ["--derivatives", "path/to/nowhere"]
    mock_basic_args: ClassVar[list[str]] = [
        "script_name",
        "path/to/input",
        "path/to/output",
        "participant",
    ]
    mock_all_args: ClassVar[list[str]] = mock_basic_args + mock_args_special

    @given(sb_st.inputs_configs())
    @allow_function_scoped
    def test_dynamic_inputs(self, mocker: MockerFixture, pybids_inputs: InputsConfig):
        p = create_parser()
        add_dynamic_args(p, copy.deepcopy(parse_args), pybids_inputs)
        magic_filters = list(
            it.chain.from_iterable(
                [[f"--filter-{key}", "entity=value"] for key in pybids_inputs]
            )
        )
        magic_wildcards = list(
            it.chain.from_iterable(
                [[f"--wildcards-{key}", "test"] for key in pybids_inputs]
            )
        )
        magic_path = list(
            it.chain.from_iterable([[f"--path-{key}", "test"] for key in pybids_inputs])
        )
        mocker.patch.object(
            sys,
            "argv",
            self.mock_all_args + magic_filters + magic_wildcards + magic_path,
        )
        args = parse_snakebids_args(p)
        for key in pybids_inputs:
            key_identifier = key.replace("-", "_")  # argparse does this
            assert isinstance(args.args_dict[f"path_{key_identifier}"], str)
            assert isinstance(args.args_dict[f"filter_{key_identifier}"], dict)
            assert isinstance(args.args_dict[f"wildcards_{key_identifier}"], list)

    @given(
        pybids_inputs=sb_st.inputs_configs(),
        flag=st.from_regex(
            re.compile(r"(?:required)|(?:any)", re.IGNORECASE), fullmatch=True
        ),
    )
    @allow_function_scoped
    def test_required_filters(
        self, mocker: MockerFixture, pybids_inputs: InputsConfig, flag: str
    ):
        p = create_parser()
        add_dynamic_args(p, copy.deepcopy(parse_args), pybids_inputs)
        magic_filters = list(
            it.chain.from_iterable(
                [[f"--filter-{key}", f"entity:{flag}"] for key in pybids_inputs]
            )
        )
        mocker.patch.object(sys, "argv", self.mock_all_args + magic_filters)

        args = parse_snakebids_args(p)
        for key in pybids_inputs:
            key_identifier = key.replace("-", "_")
            assert args.args_dict[f"filter_{key_identifier}"]["entity"] is True

    @given(
        pybids_inputs=sb_st.inputs_configs(),
        flag=st.from_regex(re.compile(r"optional", re.IGNORECASE), fullmatch=True),
    )
    @allow_function_scoped
    def test_optional_filters(
        self, mocker: MockerFixture, pybids_inputs: InputsConfig, flag: str
    ):
        p = create_parser()
        add_dynamic_args(p, copy.deepcopy(parse_args), pybids_inputs)
        magic_filters = list(
            it.chain.from_iterable(
                [[f"--filter-{key}", f"entity:{flag}"] for key in pybids_inputs]
            )
        )
        mocker.patch.object(sys, "argv", self.mock_all_args + magic_filters)

        args = parse_snakebids_args(p)
        for key in pybids_inputs:
            key_identifier = key.replace("-", "_")
            assert (
                args.args_dict[f"filter_{key_identifier}"]["entity"] is OptionalFilter
            )

    @given(
        pybids_inputs=sb_st.inputs_configs(),
        flag=st.from_regex(re.compile(r"none", re.IGNORECASE), fullmatch=True),
    )
    @allow_function_scoped
    def test_none_filters(
        self, mocker: MockerFixture, pybids_inputs: InputsConfig, flag: str
    ):
        p = create_parser()
        add_dynamic_args(p, copy.deepcopy(parse_args), pybids_inputs)
        magic_filters = list(
            it.chain.from_iterable(
                [[f"--filter-{key}", f"entity:{flag}"] for key in pybids_inputs]
            )
        )
        mocker.patch.object(sys, "argv", self.mock_all_args + magic_filters)

        args = parse_snakebids_args(p)
        for key in pybids_inputs:
            key_identifier = key.replace("-", "_")
            assert args.args_dict[f"filter_{key_identifier}"]["entity"] is False

    def test_fails_if_missing_arguments(
        self, parser: ArgumentParser, mocker: MockerFixture
    ):
        mocker.patch.object(sys, "argv", ["script_name"])
        with pytest.raises(SystemExit):
            parser.parse_args()

    def test_succeeds_if_given_positional_args(
        self, parser: ArgumentParser, mocker: MockerFixture
    ):
        mocker.patch.object(sys, "argv", self.mock_basic_args)
        assert isinstance(parser.parse_args(), Namespace)

    def test_converts_type_path_into_pathlike(
        self, parser: ArgumentParser, mocker: MockerFixture
    ):
        mocker.patch.object(sys, "argv", self.mock_all_args)
        args = parser.parse_args()
        assert isinstance(args.derivatives[0], Path)

    def test_fails_if_undefined_type_given(self):
        parse_args_copy = copy.deepcopy(parse_args)
        parse_args_copy["--new-param"] = {
            "help": "Generic Help Message",
            "type": "UnheardClass",
        }
        with pytest.raises(ConfigError, match="could not be resolved"):
            add_dynamic_args(create_parser(), parse_args_copy, pybids_inputs)

    def test_convert_arg_to_builtin(
        self, parser: ArgumentParser, mocker: MockerFixture
    ):
        new_args = {
            "--new-param": {
                "help": "Generic Help Message",
                "type": "int",
            }
        }
        mocker.patch.object(sys, "argv", [*self.mock_all_args, "--new-param", "12"])
        add_dynamic_args(parser, new_args, {})
        args = parser.parse_args()
        assert isinstance(args.new_param, int)

    def test_non_serialiable_type_raises_error(self, parser: ArgumentParser):
        new_args = {
            "--new-param": {
                "help": "Generic Help Message",
                "type": "snakebids.utils.utils.BidsEntity",
            }
        }
        with pytest.raises(ConfigError, match="cannot be serialized into yaml"):
            add_dynamic_args(parser, new_args, {})

    def test_using_module_as_type_gives_error(self, parser: ArgumentParser):
        new_args = {
            "--new-param": {
                "help": "Generic Help Message",
                "type": "snakebids.utils.utils",
            }
        }
        with pytest.raises(ConfigError, match="cannot be used as a type"):
            add_dynamic_args(parser, new_args, {})

    def test_using_class_method_as_type_gives_error(self, parser: ArgumentParser):
        new_args = {
            "--new-param": {
                "help": "Generic Help Message",
                "type": "snakebids.utils.utils.BidsEntity.from_tag",
            }
        }
        with pytest.raises(ConfigError, match="could not be resolved"):
            add_dynamic_args(parser, new_args, {})

    def test_resolves_paths(self, parser: ArgumentParser, mocker: MockerFixture):
        mocker.patch.object(sys, "argv", self.mock_all_args)
        args = parse_snakebids_args(parser).args_dict
        assert args["derivatives"][0] == Path.cwd() / "path/to/nowhere"


def test_dash_syntax_in_config_cli_args(parser: ArgumentParser, mocker: MockerFixture):
    mocker.patch.object(
        sys,
        "argv",
        [
            "script_name",
            "path/to/input",
            "path/to/output",
            "participant",
            "--participant-label",
            "12345",
            "--arg_using_dash_syntax",
            "7890",
        ],
    )
    args = parse_snakebids_args(parser)
    assert args.args_dict["participant_label"][0] == "12345"
    assert args.args_dict["arg_using_dash_syntax"][0] == "7890"
