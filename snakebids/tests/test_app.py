from __future__ import annotations

import copy
import json
import sys
from importlib import metadata
from pathlib import Path
from typing import Any, cast

import hypothesis.strategies as st
import pytest
import snakemake
from hypothesis import HealthCheck, assume, example, given, settings
from pytest_mock.plugin import MockerFixture

from snakebids.app import update_config
from snakebids.cli import SnakebidsArgs
from snakebids.tests import strategies as sb_st
from snakebids.types import InputConfig, InputsConfig, OptionalFilter

from .. import app as sn_app
from ..app import SnakeBidsApp
from .mock.config import config


@pytest.fixture
def app(mocker: MockerFixture):
    app = SnakeBidsApp(
        Path("app"),
        skip_parse_args=False,
        snakefile_path=Path("Snakefile"),
        configfile_path=Path("mock/config.yaml"),
        config=copy.deepcopy(config),
    )
    app.config["analysis_level"] = "participant"
    app.config["snakemake_args"] = []
    app.config["pybidsdb_reset"] = False
    mocker.patch.object(sn_app, "update_config", return_value=app.config)
    return app


class TestUpdateConfig:
    @given(
        input_config=sb_st.input_configs(),
        drop_wildcards=st.booleans(),
    )
    def test_magic_args(
        self,
        input_config: InputConfig,
        drop_wildcards: bool,
    ):
        config_copy: dict[str, Any] = copy.deepcopy(config)
        config_copy["bids_dir"] = "root"
        config_copy["output_dir"] = "app"
        if drop_wildcards:
            del config_copy["pybids_inputs"]["bold"]["wildcards"]
        args = SnakebidsArgs(
            force=False,
            outputdir=Path("app"),
            snakemake_args=[],
            args_dict={
                "filter_bold": input_config.get("filters"),
                "wildcards_bold": input_config.get("wildcards"),
                "path_bold": input_config.get("custom_path"),
            },
        )
        update_config(config_copy, args)
        inputs_config: InputsConfig = cast(InputsConfig, config_copy["pybids_inputs"])
        if "filters" in input_config:
            for key, value in input_config["filters"].items():
                assert inputs_config["bold"].get("filters", {}).get(key) == value
        if "wildcards" in input_config:
            assert set(input_config["wildcards"]) <= set(
                inputs_config["bold"].get("wildcards", [])
            )
        if "custom_path" in input_config:
            assert (
                inputs_config["bold"].get("custom_path", "")
                == Path(input_config["custom_path"]).resolve()
            )

    @given(
        inputs_config=sb_st.inputs_configs(),
    )
    def test_magic_optional_filter(
        self,
        inputs_config: InputsConfig,
    ):
        config_copy: dict[str, Any] = copy.deepcopy(config)
        config_copy["bids_dir"] = "root"
        config_copy["output_dir"] = "app"
        config_copy["pybids_inputs"] = inputs_config
        args_dict: dict[str, Any] = {
            f"filter_{input_}": {
                entity: OptionalFilter for entity in value.get("filters", {})
            }
            for input_, value in inputs_config.items()
        }
        for input_ in inputs_config:
            args_dict[f"wildcards_{input_}"] = None
            args_dict[f"path_{input_}"] = None
        args = SnakebidsArgs(
            force=False,
            outputdir=Path("app"),
            snakemake_args=[],
            args_dict=args_dict,
        )
        update_config(config_copy, args)
        inputs_config = config_copy["pybids_inputs"]
        for input_config in inputs_config.values():
            assert len(input_config.get("filters", {})) == 0


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
            skip_parse_args=False,
            snakefile_path=Path("Snakefile"),
            configfile_path=Path("mock/config.yaml"),
            config=copy.deepcopy(config),
        )
        app.config["analysis_level"] = "participant"
        app.config["snakemake_args"] = []
        mocker.patch.object(
            sn_app,
            "update_config",
            side_effect=lambda config, sn_args: (  # type: ignore
                config.update(sn_args.args_dict)  # type: ignore
            ),
        )

        # Prepare expected config
        expected_config = copy.deepcopy(app.config)
        expected_config.update(
            {
                "root": "",
                "snakemake_dir": Path("app").resolve(),
                "pybidsdb_dir": Path("/tmp/output/.db"),
                "pybidsdb_reset": True,
                "snakefile": Path("Snakefile"),
                "output_dir": outputdir.resolve(),
                "snakemake_version": metadata.version("snakemake"),
                "snakebids_version": metadata.version("snakebids"),
                "app_version": "unknown",  # not installing a snakebids app here
            }
        )
        if root == "app" and not tail:
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
            pybidsdb_reset=True,
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

    def test_pybidsdb_path_resolved(self, mocker: MockerFixture):
        self.io_mocks(mocker)
        mocker.patch.object(
            sys,
            "argv",
            ["./run.sh", "input", "output", "participant", "--pybidsdb-dir", ".pybids"],
        )

        # Prepare app and initial config values
        app = SnakeBidsApp(
            Path("app"),
            skip_parse_args=False,
            snakefile_path=Path("Snakefile"),
            configfile_path=Path("mock/config.yaml"),
            config=copy.deepcopy(config),
        )

        # Prepare expected config
        try:
            app.run_snakemake()
        except SystemExit as e:
            print("System exited prematurely")
            print(e)

        assert app.config["pybidsdb_dir"] == Path(".pybids").resolve()

    def test_plugin_args(self, mocker: MockerFixture, app: SnakeBidsApp):
        """Test that plugins have access to args parsed from the CLI."""
        # Get mocks for all the io functions
        self.io_mocks(mocker)
        mocker.patch.object(
            sn_app,
            "update_config",
            side_effect=(
                lambda config, sn_args: config.update(sn_args.args_dict)  # type: ignore
            ),
        )
        mocker.patch.object(
            sys, "argv", ["script_name", "path/to/input", "app/results", "participant"]
        )

        def plugin(my_app: SnakeBidsApp):
            my_app.foo = my_app.args.outputdir  # type: ignore

        app.plugins.extend([plugin])
        try:
            app.run_snakemake()
        except SystemExit as e:
            print("System exited prematurely")
            print(e)

        assert app.foo == Path("app/results").resolve()  # type: ignore

    def test_plugins(self, mocker: MockerFixture, app: SnakeBidsApp):
        # Get mocks for all the io functions
        self.io_mocks(mocker)
        mocker.patch.object(
            sn_app,
            "update_config",
            side_effect=(
                lambda config, sn_args: config.update(sn_args.args_dict)  # type: ignore
            ),
        )
        output_dir = Path("app") / "results"
        app.args = SnakebidsArgs(
            force=False,
            outputdir=output_dir,
            snakemake_args=[],
            args_dict={"output_dir": output_dir.resolve()},
        )

        def plugin(my_app: SnakeBidsApp):
            my_app.foo = "bar"  # type: ignore

        app.plugins.extend([plugin])
        try:
            app.run_snakemake()
        except SystemExit as e:
            print("System exited prematurely")
            print(e)

        assert app.foo == "bar"  # type: ignore

    def test_get_app_version_no_package(self, app: SnakeBidsApp):
        assert app.version is None

    def test_get_app_version_package(self, mocker: MockerFixture):
        metadata_pkg = (
            "importlib.metadata" if sys.version_info >= (3, 8) else "importlib_metadata"
        )
        mock = mocker.patch(f"{metadata_pkg}.version", return_value="0.1.0")
        app = SnakeBidsApp(
            Path("my_app"),
            snakefile_path=Path("Snakefile"),
            configfile_path=Path("mock/config.yaml"),
            config=copy.deepcopy(config),
        )

        assert app.version == "0.1.0"
        mock.assert_called_once_with("my_app")


class TestGenBoutiques:
    def test_boutiques_descriptor(self, tmp_path: Path, app: SnakeBidsApp):
        descriptor_path = tmp_path / "descriptor.json"
        app.create_descriptor(descriptor_path)
        with open(descriptor_path, encoding="utf-8") as descriptor_file:
            descriptor_json = json.load(descriptor_file)
            assert "command-line" in descriptor_json
            assert "inputs" in descriptor_json
