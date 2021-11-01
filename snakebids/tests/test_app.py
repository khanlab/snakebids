# Core
import sys, copy
from pathlib import Path

# Testing 
import pytest

# Typings
from argparse import ArgumentParser
from os import PathLike
from configargparse import Namespace

# Fixtures
from pytest_mock.plugin import MockerFixture

# Local
from ..app import SnakeBidsApp, resolve_path
from .mock.config import config

def init_snakebids_app(self):
    self.configfile_path="mock/config.yaml"
    self.snakefile = "mock/Snakefile"
    self.config = copy.deepcopy(config)

class TestResolvePath:
    @pytest.fixture
    def arg_dict(self):
        return {
            "bids_dir": "path/to/input",
            "output_dir": "path/to/output",
            "analysis_level": "participant",
            "--derivatives": [
                "path/to/deriv1",
                "path/to/deriv2"
            ]
        }

    def test_does_not_change_dict_without_paths(self, arg_dict):
        arg_dict_copy = copy.deepcopy(arg_dict)
        resolved = {key: resolve_path(value) for key, value in arg_dict.items()}
        assert resolved == arg_dict_copy
    
    def test_resolves_all_paths(self, arg_dict):
        arg_dict["--derivatives"] = [
            Path("path/to/deriv1"),
            Path("path/to/deriv2")
        ]
        arg_dict_copy = copy.deepcopy(arg_dict)
        arg_dict_copy["--derivatives"] = [
            p.resolve() for p in arg_dict_copy["--derivatives"]
        ]
        resolved = {key: resolve_path(value) for key, value in arg_dict.items()}
        assert resolved == arg_dict_copy

class TestArgTypeAnnotation:
    mock_args_special= ["--derivatives", "path/to/nowhere"]
    mock_basic_args = ["script_name", "path/to/input", "path/to/output", "participant"]
    mock_all_args = mock_basic_args + mock_args_special 

    @pytest.fixture
    def app(self, mocker: MockerFixture):
        mocker.patch.object(SnakeBidsApp, '__init__', init_snakebids_app)
        return SnakeBidsApp()

    @pytest.fixture
    def parser(self, app):
        return app._SnakeBidsApp__create_parser()

    def test_snakebids_app_is_properly_mocked(self, app):
        assert isinstance(app, SnakeBidsApp)
        assert not hasattr(app, "parser")

    def test_fails_if_missing_arguments(self, parser: ArgumentParser, mocker: MockerFixture):
        mocker.patch.object(sys, 'argv', ["script_name"])
        with pytest.raises(SystemExit):
            parser.parse_args()

    def test_succeeds_if_given_positional_args(self, parser: ArgumentParser, mocker: MockerFixture):
        mocker.patch.object(sys, 'argv', self.mock_basic_args)
        assert isinstance(parser.parse_args(), Namespace)

    def test_converts_type_path_into_pathlike(self, parser: ArgumentParser, mocker: MockerFixture):
        mocker.patch.object(sys, 'argv', self.mock_all_args)
        args = parser.parse_args()
        assert isinstance(getattr(args, "derivatives")[0], PathLike)
    
    def test_fails_if_undefined_type_given(self, app: SnakeBidsApp, mocker: MockerFixture):
        app.config["parse_args"]["--new-param"] = {
            "help": "Generic Help Message",
            "type": "UnheardClass"
        }
        with pytest.raises(TypeError):
            app._SnakeBidsApp__create_parser()

    def test_resolves_paths(self, app: SnakeBidsApp, mocker: MockerFixture):
        mocker.patch.object(sys, 'argv', self.mock_all_args)
        app.parser = app._SnakeBidsApp__create_parser()
        app._SnakeBidsApp__parse_args()
        assert app.config["bids_dir"] == Path.cwd() / "path/to/input"
        assert app.config["derivatives"][0] == Path.cwd() / "path/to/nowhere" 
