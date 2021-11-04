# Core
import sys, copy
from pathlib import Path
from typing import Dict
from unittest.mock import MagicMock

# Testing 
import pytest

# Typings
from argparse import ArgumentParser
from os import PathLike
from configargparse import Namespace

import snakemake
# Fixtures
from pytest_mock.plugin import MockerFixture

# Local
from .. import app as sn_app
from ..app import SnakeBidsApp, resolve_path

from .mock.config import config


def init_snakebids_app(self):
    self.configfile_path="mock/config.yaml"
    self.snakefile = "mock/Snakefile"
    self.retrofit = False
    self.config = copy.deepcopy(config)
    self.config["analysis_level"] = "participant"
    self.config["snakemake_args"] = []

@pytest.fixture
def app(mocker: MockerFixture):
    mocker.patch.object(SnakeBidsApp, '__init__', init_snakebids_app)
    return SnakeBidsApp()

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
    def parser(self, app):
        return app._create_parser()

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
            app._create_parser()

    def test_resolves_paths(self, app: SnakeBidsApp, mocker: MockerFixture):
        mocker.patch.object(sys, 'argv', self.mock_all_args)
        app.parser = app._create_parser()
        app._parse_args()
        assert app.config["bids_dir"] == Path.cwd() / "path/to/input"
        assert app.config["derivatives"][0] == Path.cwd() / "path/to/nowhere" 


class TestRunSnakemake:
    @pytest.fixture
    def io_mocks(self, mocker: MockerFixture):
        return {
            "write_output_mode": mocker.patch.object(sn_app, 'write_output_mode'),
            "prepare_output": mocker.patch.object(sn_app, 'prepare_output'),
            "write_config": mocker.patch.object(sn_app, 'write_config_file'),
            "snakemake": mocker.patch.object(snakemake, 'main'),
        }

    def test_runs_in_workflow_mode(
        self, io_mocks: Dict[str, MagicMock], app: SnakeBidsApp
    ):
        app.workflow_mode = True
        app.outputdir = Path("/tmp/output")
        app.snakemake_dir = Path("app")
        app.force = False
        expected_config = copy.deepcopy(app.config)
        expected_config["output_dir"] = "/tmp/output/results"

        io_mocks["prepare_output"].return_value = Path("/tmp/output/results")

        try:
            app.run_snakemake()
        except SystemExit as e:
            print("System exited prematurely")
            print(e)
        
        io_mocks["write_output_mode"].assert_not_called()
        io_mocks["prepare_output"].assert_called_once_with(
            Path("app"),
            Path("/tmp/output"),
            "workflow",
            False
        )
        io_mocks["write_config"].assert_called_once_with(
            Path("/tmp/output/mock/config.yaml"),
            expected_config,
            False
        )
        io_mocks["snakemake"].assert_called_once_with([
            "--snakefile",
            app.snakefile,
            "--directory",
            "/tmp/output",
            "--configfile",
            "/tmp/output/mock/config.yaml"
        ])

    def test_runs_in_bidsapp_mode(
        self, io_mocks: Dict[str, MagicMock], app: SnakeBidsApp
    ):
        app.workflow_mode = False
        app.outputdir = Path("/tmp/output")
        app.snakemake_dir = Path("app")
        app.force = False
        expected_config = copy.deepcopy(app.config)

        io_mocks["prepare_output"].return_value = Path("/tmp/output")

        try:
            app.run_snakemake()
        except SystemExit as e:
            print("System exited prematurely")
            print(e)
        
        io_mocks["write_output_mode"].assert_not_called()
        io_mocks["prepare_output"].assert_called_once_with(
            Path("app"),
            Path("/tmp/output"),
            "bidsapp",
            False
        )
        io_mocks["write_config"].assert_called_once_with(
            Path("/tmp/output/code/config.yaml"),
            expected_config,
            True
        )
        io_mocks["snakemake"].assert_called_once_with([
            "--snakefile",
            app.snakefile,
            "--directory",
            "app",
            "--configfile",
            "/tmp/output/code/config.yaml"
        ])

    def test_runs_in_workflow_mode_when_output_same_as_snakebids_app(
        self,
        io_mocks: Dict[str, MagicMock],
        app: SnakeBidsApp
    ):
        app.workflow_mode = False
        app.outputdir = Path("app")
        app.snakemake_dir = Path("app")
        app.force = False
        expected_config = copy.deepcopy(app.config)
        expected_config["output_dir"] = "app/results"

        io_mocks["prepare_output"].return_value = Path("app/results")

        try:
            app.run_snakemake()
        except SystemExit as e:
            print("System exited prematurely")
            print(e)
        
        io_mocks["write_output_mode"].assert_called_once_with(
            Path("app/.snakebids"),
            "workflow"
        )
        io_mocks["prepare_output"].assert_called_once_with(
            Path("app"),
            Path("app"),
            "workflow",
            False
        )
        io_mocks["write_config"].assert_called_once_with(
            Path("app/mock/config.yaml"),
            expected_config,
            True
        )
        io_mocks["snakemake"].assert_called_once_with([
            "--snakefile",
            app.snakefile,
            "--directory",
            "app",
            "--configfile",
            str(Path("app/mock/config.yaml").resolve())
        ])
