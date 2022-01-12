# pylint: disable=protected-access, redefined-outer-name
from __future__ import absolute_import

import copy
import sys
from argparse import ArgumentParser
from os import PathLike
from pathlib import Path

import pytest
from configargparse import Namespace
from pytest_mock.plugin import MockerFixture

from snakebids.cli import (
    _resolve_path,
    add_dynamic_args,
    create_parser,
    parse_snakebids_args,
)

from .mock.config import parse_args, pybids_inputs


@pytest.fixture
def parser():
    p, p_run = create_parser()
    add_dynamic_args(p_run, copy.deepcopy(parse_args), copy.deepcopy(pybids_inputs))
    return p


class TestResolvePath:
    @pytest.fixture
    def arg_dict(self):
        return {
            "bids_dir": "path/to/input",
            "output_dir": "path/to/output",
            "analysis_level": "participant",
            "--derivatives": ["path/to/deriv1", "path/to/deriv2"],
        }

    def test_does_not_change_dict_without_paths(self, arg_dict):
        arg_dict_copy = copy.deepcopy(arg_dict)
        resolved = {key: _resolve_path(value) for key, value in arg_dict.items()}
        assert resolved == arg_dict_copy

    def test_resolves_all_paths(self, arg_dict):
        arg_dict["--derivatives"] = [Path("path/to/deriv1"), Path("path/to/deriv2")]
        arg_dict_copy = copy.deepcopy(arg_dict)
        arg_dict_copy["--derivatives"] = [
            p.resolve() for p in arg_dict_copy["--derivatives"]
        ]
        resolved = {key: _resolve_path(value) for key, value in arg_dict.items()}
        assert resolved == arg_dict_copy


class TestAddDynamicArgs:
    mock_args_special = ["--derivatives", "path/to/nowhere"]
    mock_basic_args = [
        "script_name",
        "run",
        "path/to/input",
        "path/to/output",
        "participant",
    ]
    mock_all_args = mock_basic_args + mock_args_special

    def test_fails_if_missing_arguments(
        self, parser: ArgumentParser, mocker: MockerFixture
    ):
        mocker.patch.object(sys, "argv", ["run", "script_name"])
        with pytest.raises(SystemExit):
            parser.parse_args()

    def test_succeeds_if_given_positional_args(
        self, parser: ArgumentParser, mocker: MockerFixture
    ):
        mocker.patch.object(sys, "argv", self.mock_basic_args)
        assert isinstance(parser.parse_args(), Namespace)

    def test_succeeds_if_given_boutiques_args(
        self, parser: ArgumentParser, mocker: MockerFixture
    ):
        mocker.patch.object(sys, "argv", ["script_name", "boutiques", "descriptor.txt"])
        assert isinstance(parser.parse_args(), Namespace)

    def test_converts_type_path_into_pathlike(
        self, parser: ArgumentParser, mocker: MockerFixture
    ):
        mocker.patch.object(sys, "argv", self.mock_all_args)
        args = parser.parse_args()
        assert isinstance(getattr(args, "derivatives")[0], PathLike)

    def test_fails_if_undefined_type_given(self):
        parse_args_copy = copy.deepcopy(parse_args)
        parse_args_copy["--new-param"] = {
            "help": "Generic Help Message",
            "type": "UnheardClass",
        }
        with pytest.raises(TypeError):
            add_dynamic_args(create_parser()[1], parse_args_copy, pybids_inputs)

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
            "run",
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
