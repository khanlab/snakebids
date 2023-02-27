# pylint: disable=protected-access, redefined-outer-name
from __future__ import absolute_import

import copy
import json
from pathlib import Path

import hypothesis.strategies as st
import pytest
import snakemake
from hypothesis import HealthCheck, assume, example, given, settings
from pytest_mock.plugin import MockerFixture

from snakebids.cli import SnakebidsArgs

from .. import app as sn_app
from ..app import SnakeBidsApp
from .mock.config import config


@pytest.fixture
def app(mocker: MockerFixture):
    app = SnakeBidsApp(
        Path("app"),
        False,
        snakefile_path=Path("Snakefile"),
        configfile_path=Path("mock/config.yaml"),
        config=copy.deepcopy(config),
    )
    app.config["analysis_level"] = "participant"
    app.config["snakemake_args"] = []
    app.config["pybids_db_reset"] = False
    mocker.patch.object(sn_app, "update_config", return_value=app.config)
    return app


class TestRunSnakemake:
    valid_chars = st.characters(
        min_codepoint=48, max_codepoint=122, whitelist_categories=["Ll", "Lu"]
    )

    def io_mocks(self, mocker: MockerFixture):
        return {
            "write_output_mode": mocker.patch.object(sn_app, "write_output_mode"),
            "prepare_output": mocker.patch.object(sn_app, "prepare_bidsapp_output"),
            "write_config": mocker.patch.object(sn_app, "write_config_file"),
            "snakemake": mocker.patch.object(snakemake, "main"),
        }

    @given(
        root=st.one_of(st.sampled_from(["app", "app/results"]), valid_chars),
        tail=valid_chars,
    )
    @example(root="app", tail="")
    @example(root="app/results", tail="")
    # mocker is function_scoped, but doesn't cause any problems here
    @settings(suppress_health_check=[HealthCheck.function_scoped_fixture])
    def test_runs_in_correct_mode(
        self,
        mocker: MockerFixture,
        root: str,
        tail: str,
    ):
        # Test any valid path. Whether the paths are bytes is irrelevant, and only
        # causes problems
        assume(not isinstance(root, bytes))
        assume(not isinstance(tail, bytes))
        # Prevent the case where tail == "."
        assume(Path(root).resolve() / tail != Path(root).resolve())

        # Get mocks for all the io functions
        io_mocks = self.io_mocks(mocker)

        # Compose outputdir from root and tail
        outputdir = Path(root) / tail

        # Prepare app and initial config values
        app = SnakeBidsApp(
            Path("app"),
            False,
            snakefile_path=Path("Snakefile"),
            configfile_path=Path("mock/config.yaml"),
            config=copy.deepcopy(config),
        )
        app.config["analysis_level"] = "participant"
        app.config["snakemake_args"] = []
        mocker.patch.object(
            sn_app,
            "update_config",
            side_effect=lambda config, sn_args: config.update(sn_args.args_dict),
        )

        # Prepare expected config
        expected_config = copy.deepcopy(app.config)
        expected_config.update(
            {
                "root": "",
                "snakemake_dir": Path("app").resolve(),
                "pybids_db_dir": Path("/tmp/output/.db"),
                "pybids_db_reset": True,
                "snakefile": Path("Snakefile"),
                "output_dir": outputdir.resolve(),
                "skip_bids_validation": False,
            }
        )
        if root == "app" and tail == "":
            expected_config["output_dir"] /= "results"
            expected_config["root"] = "results"

        app.args = SnakebidsArgs(
            force=False,
            outputdir=outputdir,
            snakemake_args=[],
            # This Dict is necessary for updating config, since "update_config" is
            # patched
            args_dict={"output_dir": outputdir.resolve()},
            pybidsdb_dir=Path("/tmp/output/.db"),
            reset_db=True,
        )

        try:
            app.run_snakemake()
        except SystemExit as e:
            print("System exited prematurely")
            print(e)

        # First condition: outputdir is an arbitrary path
        if root not in ["app", "app/results"] or (root == "app" and tail):
            cwd = outputdir.resolve()
            new_config = outputdir.resolve() / "mock/config.yaml"
            io_mocks["write_output_mode"].assert_not_called()
            io_mocks["prepare_output"].assert_called_once_with(
                outputdir.resolve(), False
            )
            io_mocks["write_config"].assert_called_once_with(
                config_file=new_config,
                data=expected_config,
                force_overwrite=True,
            )

        # Second condition: outputdir is equal to snakemake_dir, or under
        # snakemake_dir/results
        else:
            new_config = Path("app").resolve() / "mock/config.yaml"
            cwd = Path("app").resolve()
            io_mocks["write_output_mode"].assert_called_once_with(
                Path("app/.snakebids").resolve(), "workflow"
            )
            io_mocks["write_config"].assert_called_once_with(
                config_file=new_config,
                data=expected_config,
                force_overwrite=True,
            )

        io_mocks["snakemake"].assert_called_once_with(
            [
                "--snakefile",
                str(app.snakefile_path),
                "--directory",
                str(cwd),
                "--configfile",
                str(new_config),
            ]
        )

    def test_plugins(self, mocker: MockerFixture, app: SnakeBidsApp):
        # Get mocks for all the io functions
        self.io_mocks(mocker)
        mocker.patch.object(
            sn_app,
            "update_config",
            side_effect=lambda config, sn_args: config.update(sn_args.args_dict),
        )
        output_dir = Path("app") / "results"
        app.args = SnakebidsArgs(
            force=False,
            outputdir=output_dir,
            snakemake_args=[],
            args_dict={"output_dir": output_dir.resolve()},
        )

        def plugin(my_app):
            my_app.foo = "bar"

        app.add_plugins([plugin])
        try:
            app.run_snakemake()
        except SystemExit as e:
            print("System exited prematurely")
            print(e)

        assert app.foo == "bar"


class TestGenBoutiques:
    def test_boutiques_descriptor(self, tmp_path: Path, app: SnakeBidsApp):
        descriptor_path = tmp_path / "descriptor.json"
        app.create_descriptor(descriptor_path)
        with open(descriptor_path, encoding="utf-8") as descriptor_file:
            descriptor_json = json.load(descriptor_file)
            assert "command-line" in descriptor_json
            assert "inputs" in descriptor_json
