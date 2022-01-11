# pylint: disable=protected-access, redefined-outer-name
from __future__ import absolute_import

import copy
from pathlib import Path
from typing import Dict
from unittest.mock import MagicMock

import pytest
import snakemake
from pytest_mock.plugin import MockerFixture

from snakebids.cli import SnakebidsArgs

from .. import app as sn_app
from ..app import SnakeBidsApp
from ..cli import create_parser
from .mock.config import config


@pytest.fixture
def app(mocker: MockerFixture):
    parser, parser_run = create_parser()
    app = SnakeBidsApp(
        snakemake_dir=Path("app"),
        skip_parse_args=False,
        parser=parser,
        parser_run=parser_run,
        snakefile_path=Path("Snakefile"),
        configfile_path=Path("mock/config.yaml"),
        config=copy.deepcopy(config),
    )
    app.config["analysis_level"] = "participant"
    app.config["snakemake_args"] = []
    mocker.patch.object(sn_app, "update_config", return_value=app.config)
    return app


def test_update_config(app: SnakeBidsApp, mocker: MockerFixture):
    pass


class TestRunSnakemake:
    @pytest.fixture
    def io_mocks(self, mocker: MockerFixture):
        return {
            "write_output_mode": mocker.patch.object(sn_app, "write_output_mode"),
            "prepare_output": mocker.patch.object(sn_app, "prepare_output"),
            "write_config": mocker.patch.object(sn_app, "write_config_file"),
            "snakemake": mocker.patch.object(snakemake, "main"),
        }

    def test_runs_in_workflow_mode(
        self, io_mocks: Dict[str, MagicMock], app: SnakeBidsApp
    ):
        expected_config = copy.deepcopy(app.config)
        expected_config["output_dir"] = "/tmp/output/results"
        expected_config["root"] = "results"
        expected_config["snakemake_dir"] = Path("app").resolve()
        expected_config["snakefile"] = Path("Snakefile")

        io_mocks["prepare_output"].return_value = Path("/tmp/output/results")

        app.args = SnakebidsArgs(
            workflow_mode=True,
            force=False,
            outputdir=Path("/tmp/output"),
            retrofit=False,
            snakemake_args=[],
            args_dict={},
        )

        try:
            app.run_snakemake()
        except SystemExit as e:
            print("System exited prematurely")
            print(e)

        io_mocks["write_output_mode"].assert_not_called()
        io_mocks["prepare_output"].assert_called_once_with(
            Path("app").resolve(), Path("/tmp/output"), "workflow", False
        )
        io_mocks["write_config"].assert_called_once_with(
            config_file=Path("/tmp/output/mock/config.yaml"),
            data=expected_config,
            force_overwrite=False,
        )
        io_mocks["snakemake"].assert_called_once_with(
            [
                "--snakefile",
                str(app.snakefile_path),
                "--directory",
                "/tmp/output",
                "--configfile",
                "/tmp/output/mock/config.yaml",
            ]
        )

    def test_runs_in_bidsapp_mode(
        self, io_mocks: Dict[str, MagicMock], app: SnakeBidsApp
    ):
        expected_config = copy.deepcopy(app.config)
        expected_config["root"] = ""
        expected_config["snakemake_dir"] = Path("app").resolve()
        expected_config["snakefile"] = Path("Snakefile")

        io_mocks["prepare_output"].return_value = Path("/tmp/output")

        app.args = SnakebidsArgs(
            workflow_mode=False,
            force=False,
            outputdir=Path("/tmp/output"),
            retrofit=False,
            snakemake_args=[],
            args_dict={},
        )

        try:
            app.run_snakemake()
        except SystemExit as e:
            print("System exited prematurely")
            print(e)

        io_mocks["write_output_mode"].assert_not_called()
        io_mocks["prepare_output"].assert_called_once_with(
            Path("app").resolve(), Path("/tmp/output"), "bidsapp", False
        )
        io_mocks["write_config"].assert_called_once_with(
            config_file=Path("/tmp/output/mock/config.yaml"),
            data=expected_config,
            force_overwrite=True,
        )
        io_mocks["snakemake"].assert_called_once_with(
            [
                "--snakefile",
                str(app.snakefile_path),
                "--directory",
                "/tmp/output",
                "--configfile",
                "/tmp/output/mock/config.yaml",
            ]
        )

    def test_runs_in_workflow_mode_when_output_same_as_snakebids_app(
        self, io_mocks: Dict[str, MagicMock], app: SnakeBidsApp
    ):
        expected_config = copy.deepcopy(app.config)
        expected_config["output_dir"] = "app/results"
        expected_config["root"] = "results"
        expected_config["snakemake_dir"] = Path("app").resolve()
        expected_config["snakefile"] = Path("Snakefile")

        io_mocks["prepare_output"].return_value = Path("app/results")

        app.args = SnakebidsArgs(
            workflow_mode=False,
            force=False,
            outputdir=Path("app").resolve(),
            retrofit=False,
            snakemake_args=[],
            args_dict={},
        )

        try:
            app.run_snakemake()
        except SystemExit as e:
            print("System exited prematurely")
            print(e)

        io_mocks["write_output_mode"].assert_called_once_with(
            Path("app/.snakebids").resolve(), "workflow"
        )
        io_mocks["prepare_output"].assert_called_once_with(
            Path("app").resolve(), Path("app").resolve(), "workflow", False
        )
        io_mocks["write_config"].assert_called_once_with(
            config_file=Path("app/mock/config.yaml").resolve(),
            data=expected_config,
            force_overwrite=True,
        )
        io_mocks["snakemake"].assert_called_once_with(
            [
                "--snakefile",
                str(app.snakefile_path),
                "--directory",
                str(Path("app").resolve()),
                "--configfile",
                str(Path("app/mock/config.yaml").resolve()),
            ]
        )
