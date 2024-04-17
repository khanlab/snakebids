from __future__ import annotations

import argparse
import copy
import tempfile
from pathlib import Path
from typing import Any, Mapping, Sequence

import pytest
from hypothesis import given
from hypothesis import strategies as st
from pytest_mock import MockerFixture
from pytest_mock.plugin import MockType

from snakebids import bidsapp
from snakebids.exceptions import ConfigError, RunError
from snakebids.plugins import Version
from snakebids.plugins.snakemake import (
    CONFIGFILE_CHOICES,
    SNAKEFILE_CHOICES,
    SnakemakeBidsApp,
    _resolve_path,
)
from snakebids.tests import strategies as sb_st
from snakebids.tests.helpers import allow_function_scoped
from snakebids.utils.utils import DEPRECATION_FLAG


class TestFindSnakefileConfig:
    @allow_function_scoped
    @given(snakefile=st.sampled_from(SNAKEFILE_CHOICES))
    def test_finds_snakefile(self, fakefs_tmpdir: Path, snakefile: str):
        tmp = Path(tempfile.mkdtemp(dir=fakefs_tmpdir))
        snakefile_path = tmp / snakefile
        snakefile_path.parent.mkdir(parents=True, exist_ok=True)
        snakefile_path.touch()
        app = SnakemakeBidsApp(snakemake_dir=tmp, configfile_path=None)
        assert app.snakefile_path == snakefile_path

    def test_errors_if_no_snakefile_found(self, fakefs_tmpdir: Path):
        tmp = Path(tempfile.mkdtemp(dir=fakefs_tmpdir))
        with pytest.raises(ConfigError, match="Error: no Snakefile"):
            SnakemakeBidsApp(snakemake_dir=tmp, configfile_path=None)

    @allow_function_scoped
    @given(configfile=st.sampled_from(CONFIGFILE_CHOICES))
    def test_finds_config_file(self, fakefs_tmpdir: Path, configfile: str):
        tmp = Path(tempfile.mkdtemp(dir=fakefs_tmpdir))
        configfile_path = tmp / configfile
        configfile_path.parent.mkdir(parents=True, exist_ok=True)
        configfile_path.touch()
        app = SnakemakeBidsApp(snakemake_dir=tmp, snakefile_path=Path())
        assert app.configfile_path == configfile_path

    def test_errors_if_no_configfile_found(self, fakefs_tmpdir: Path):
        tmp = Path(tempfile.mkdtemp(dir=fakefs_tmpdir))
        with pytest.raises(ConfigError, match="Error: no config file"):
            SnakemakeBidsApp(snakemake_dir=tmp, snakefile_path=Path())


class TestAddCliArguments:
    def test_snakemake_help_arg_added(self, mocker: MockerFixture):
        import snakebids.plugins.snakemake

        mock = mocker.patch.object(snakebids.plugins.snakemake, "snakemake_main")
        smk = SnakemakeBidsApp.create_empty()
        parser = argparse.ArgumentParser()
        smk.add_cli_arguments(parser)
        parser.parse_args(["--help-snakemake"])
        mock.assert_called_once_with(["-h"])


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


class TestUpdateNamespace:
    @given(
        config=st.dictionaries(
            st.text(),
            st.text(
                st.characters(
                    blacklist_categories=["Cs"], blacklist_characters=["\x00"]
                )
            ).map(Path),
        )
    )
    def test_resolves_paths(self, config: dict[str, Path]):
        smk = SnakemakeBidsApp.create_empty()
        namespace = {
            **config,
            "output_dir": "foo",
            "bids_dir": "foo",
            "some_key": "foo",
        }
        smk.update_cli_namespace(namespace, {})
        for key, path in namespace.items():
            if key == "some_key":
                assert path == "foo"
            else:
                assert isinstance(path, Path)
                assert path.is_absolute()

    @given(targets=st.dictionaries(st.text(), st.text(), min_size=1), data=st.data())
    def test_target_gets_set(self, targets: dict[str, str], data: st.DataObject):
        smk = SnakemakeBidsApp.create_empty()
        level = data.draw(st.sampled_from(list(targets)))
        namespace = {"analysis_level": level}
        config = {"targets_by_analysis_level": targets}
        smk.update_cli_namespace(namespace, config)
        assert config["snakemake_target"] == targets[level]

    @given(value=st.sampled_from([True, False]))
    def test_force_output_is_set(self, value: bool):
        smk = SnakemakeBidsApp.create_empty()
        namespace = {"force_output": value}
        smk.update_cli_namespace(namespace, {})
        assert smk.force_output == value
        assert "force_output" not in namespace


def get_io_mocks(mocker: MockerFixture):
    import snakebids.snakemake_compat  # noqa: I001
    import snakebids.plugins.snakemake as sn_app

    mocker.stopall()
    mocker.patch.object(sn_app.impm, "version", return_value="version")
    mocker.patch.object(
        snakebids.snakemake_compat.configfile,
        "open",
        mocker.mock_open(read_data='{"a_key": "a value"}'),
    )
    return {
        "write_output_mode": mocker.patch.object(sn_app, "write_output_mode"),
        "prepare_output": mocker.patch.object(sn_app, "prepare_bidsapp_output"),
        "write_config": mocker.patch.object(sn_app, "write_config"),
        "snakemake": mocker.patch.object(sn_app, "snakemake_main"),
    }


class TestFinalizeConfig:
    def check_write_config(
        self,
        mock: MockType,
        config_file: Path,
        config: dict[str, Any],
        root: Path,
    ):
        mock.assert_called_once_with(
            config_file=config_file,
            data=dict(
                config,
                snakemake_version="version",
                snakebids_version="version",
                root=root,
                app_version="unknown",
                snakemake_dir=Path().resolve(),
                snakefile=Path("Snakefile"),
            ),
            force_overwrite=True,
        )

    @given(
        path=sb_st.paths(resolve=True).filter(
            lambda p: p != Path().resolve()
            and not str(p).startswith(str(Path("results").resolve()))
        )
    )
    @allow_function_scoped
    def test_normally_preserves_paths(
        self,
        mocker: MockerFixture,
        path: Path,
    ):
        smk = SnakemakeBidsApp.create_empty()
        io_mocks = get_io_mocks(mocker)

        config = {"output_dir": path}

        smk.finalize_config(config)

        assert config["output_dir"] == path
        assert smk.cwd == path
        io_mocks["write_output_mode"].assert_not_called()
        io_mocks["prepare_output"].assert_called_once_with(path, False)

        new_config = path.resolve() / "snakebids.yaml"
        self.check_write_config(
            io_mocks["write_config"], config_file=new_config, config=config, root=Path()
        )

    @pytest.mark.parametrize("path", [Path().resolve(), Path("results").resolve()])
    def test_output_in_app_triggers_workflow_mode(
        self,
        mocker: MockerFixture,
        path: Path,
    ):
        smk = SnakemakeBidsApp.create_empty()
        io_mocks = get_io_mocks(mocker)

        config = {"output_dir": path}

        smk.finalize_config(config)

        assert config["output_dir"] == Path("results").resolve()
        assert smk.cwd == Path().resolve()
        io_mocks["write_output_mode"].assert_called_once_with(
            Path(".snakebids").resolve(), "workflow"
        )
        io_mocks["prepare_output"].assert_not_called()

        new_config = Path("snakebids.yaml").resolve()
        self.check_write_config(
            io_mocks["write_config"],
            config_file=new_config,
            config=config,
            root=Path("results"),
        )

    @given(tail=sb_st.paths(absolute=False))
    @allow_function_scoped
    def test_output_under_results_triggers_workflow_mode(
        self,
        mocker: MockerFixture,
        tail: Path,
    ):
        smk = SnakemakeBidsApp.create_empty()
        io_mocks = get_io_mocks(mocker)
        path = Path("results", tail).resolve()

        config = {"output_dir": path}

        smk.finalize_config(config)

        assert config["output_dir"] == path
        assert smk.cwd == Path().resolve()
        io_mocks["write_output_mode"].assert_called_once_with(
            Path(".snakebids").resolve(), "workflow"
        )
        io_mocks["prepare_output"].assert_not_called()

        new_config = Path("snakebids.yaml").resolve()
        self.check_write_config(
            io_mocks["write_config"],
            config_file=new_config,
            config=config,
            root=Path("results", tail),
        )

    @given(
        output=st.sampled_from(["", "results"]).map(Path) | sb_st.paths(resolve=True),
        configfile=sb_st.paths(absolute=False),
    )
    @allow_function_scoped
    def test_relative_configfile_gets_transferred(
        self, mocker: MockerFixture, output: Path, configfile: Path
    ):
        smk = SnakemakeBidsApp.create_empty()
        smk.configfile_path = configfile
        get_io_mocks(mocker)
        config = {"output_dir": output}

        smk.finalize_config(config)

        assert smk.configfile_outpath == smk.cwd / configfile

    @given(
        output=st.sampled_from(["", "results"]).map(Path) | sb_st.paths(absolute=False),
        configfile=sb_st.paths(absolute=True),
    )
    @allow_function_scoped
    def test_non_relative_configfile_is_preserved(
        self, mocker: MockerFixture, output: Path, configfile: Path
    ):
        smk = SnakemakeBidsApp.create_empty()
        smk.configfile_path = configfile
        get_io_mocks(mocker)
        config = {"output_dir": output}

        smk.finalize_config(config)

        assert smk.configfile_outpath == configfile

    @given(
        output=st.sampled_from(["", "results"]).map(Path) | sb_st.paths(resolve=True),
        configfile=sb_st.paths(),
        config_output=sb_st.paths(),
    )
    @allow_function_scoped
    def test_specified_configfile_output_always_preserved(
        self, mocker: MockerFixture, output: Path, configfile: Path, config_output: Path
    ):
        smk = SnakemakeBidsApp.create_empty()
        smk.configfile_path = configfile
        smk.configfile_outpath = config_output
        get_io_mocks(mocker)
        config = {"output_dir": output}

        smk.finalize_config(config)

        assert smk.configfile_outpath == config_output

    @given(
        err=st.text(), output=sb_st.paths(absolute=True, min_segments=1, max_segments=1)
    )
    @allow_function_scoped
    def test_exits_on_run_error(
        self,
        mocker: MockerFixture,
        capsys: pytest.CaptureFixture[str],
        err: str,
        output: Path,
    ):
        smk = SnakemakeBidsApp.create_empty()
        mocks = get_io_mocks(mocker)
        mocks["prepare_output"].side_effect = RunError(err)
        config = {"output_dir": output.resolve()}

        with pytest.raises(SystemExit):
            smk.finalize_config(config)
        cap = capsys.readouterr()
        assert err == cap.err[:-1]


class TestVersion:
    def test_get_app_version_no_package(
        self, mocker: MockerFixture, caplog: pytest.LogCaptureFixture
    ):
        plugin = SnakemakeBidsApp(
            snakemake_dir=Path(),
            configfile_path=None,
            snakefile_path=Path("Snakefile"),
        )
        assert caplog.text == ""
        get_io_mocks(mocker)
        plugin.finalize_config(
            {"output_dir": Path("foo"), "plugins.version.version": "unknown"}
        )
        assert "App version could not be identified" in caplog.text

    def test_missing_version_changed_to_unknown(self, mocker: MockerFixture):
        plugin = SnakemakeBidsApp(
            snakemake_dir=Path(),
            configfile_path=None,
            snakefile_path=Path("Snakefile"),
        )
        mocks = get_io_mocks(mocker)
        plugin.finalize_config({"output_dir": Path("foo")})

        mocks["write_config"].assert_called_once_with(
            config_file=Path("foo/snakebids.yaml"),
            data={
                "output_dir": Path("foo"),
                "snakemake_version": "version",
                "snakebids_version": "version",
                "root": Path(),
                "app_version": "unknown",
                "snakemake_dir": Path().resolve(),
                "snakefile": Path("Snakefile"),
            },
            force_overwrite=True,
        )


@pytest.mark.parametrize("root", ["app", "app/results", "foo/bar/lo", "/absolute/path"])
def test_integration(
    mocker: MockerFixture,
    root: str,
):
    # Get mocks for all the io functions
    io_mocks = get_io_mocks(mocker)

    # Prepare app and initial config values
    plugin = SnakemakeBidsApp(
        Path("app"),
        snakefile_path=Path("app/Snakefile"),
        configfile_path=Path("app/config/config.yaml"),
    )
    # dummy db_path: all io functions are mocked, so this can be arbitrary str
    db_path = "/path/to/db"
    app = bidsapp.app(
        plugins=[plugin, Version(distribution="snakebidstesting")],
    )
    app.config["targets_by_analysis_level"] = {"participant": "all"}
    app.config["snakemake_args"] = []
    app.config["pybids_inputs"] = {
        "bold": {
            "filters": {"suffix": "bold", "extension": ".nii.gz", "datatype": "func"},
            "wildcards": ["subject", "session", "acquisition", "task", "run"],
        }
    }
    app.config["parse_args"] = {
        "--additional-arg": {"nargs": "+"},
        "--derivatives": {"help": "overwrite this arg", "nargs": "*"},
    }

    # Prepare expected config
    expected_config = copy.deepcopy(app.config)
    expected_config.update(
        {
            "root": "",
            "a_key": "a value",
            "additional_arg": None,
            "snakemake_args": ["--unknown-arg"],
            "snakemake_dir": Path("app").resolve(),
            "pybidsdb_dir": Path(db_path),
            "pybidsdb_reset": True,
            "pybids_db_dir": f"{DEPRECATION_FLAG}{db_path}{DEPRECATION_FLAG}",
            "pybids_db_reset": f"{DEPRECATION_FLAG}1{DEPRECATION_FLAG}",
            "snakefile": Path("app/Snakefile"),
            "bids_dir": Path("input_dataset").resolve(),
            "output_dir": Path(root).resolve(),
            "participant_label": None,
            "exclude_participant_label": None,
            "derivatives": [],
            "snakemake_version": "version",
            "snakebids_version": "version",
            "app_version": "version",
            "plugins.version.version": "version",
            "snakemake_target": "all",
            "analysis_level": "participant",
        }
    )
    if root == "app":
        expected_config["output_dir"] /= "results"

    app.run(
        args=[
            "input_dataset",
            root,
            "participant",
            "--derivatives",
            "--pybidsdb-dir",
            db_path,
            "--pybidsdb-reset",
            "--unknown-arg",
        ],
    )

    # First condition: outputdir is equal to snakemake_dir, or under
    # snakemake_dir/results
    if root in {"app", "app/results"}:
        cwd = Path("app").resolve()
        new_config = cwd / "config/config.yaml"
        expected_config["root"] = Path("results")
        io_mocks["write_output_mode"].assert_called_once_with(
            cwd / ".snakebids", "workflow"
        )
        io_mocks["write_config"].assert_called_once_with(
            config_file=new_config,
            data=expected_config,
            force_overwrite=True,
        )

    # Second condition: outputdir is an arbitrary path
    else:
        cwd = Path(root).resolve()
        expected_config["root"] = Path()
        new_config = cwd / "config/config.yaml"
        io_mocks["write_output_mode"].assert_not_called()
        io_mocks["prepare_output"].assert_called_once_with(cwd, False)
        io_mocks["write_config"].assert_called_once_with(
            config_file=new_config,
            data=expected_config,
            force_overwrite=True,
        )

    io_mocks["snakemake"].assert_called_once_with(
        [
            "all",
            "--snakefile",
            str(plugin.snakefile_path),
            "--directory",
            str(cwd),
            "--configfile",
            str(new_config),
            "--unknown-arg",
        ]
    )
