from __future__ import absolute_import

import itertools as it
import json
import shutil
from pathlib import Path
from typing import Callable

import pytest
from pytest_mock.plugin import MockerFixture

import snakebids.output as output
from snakebids.exceptions import RunError

dirlen: Callable[[Path], int] = lambda f: len([*f.iterdir()])


@pytest.fixture
def fake_snakemake(tmp_path: Path):
    app1 = tmp_path / "app1"
    app1.mkdir()
    (app1 / "random_file").touch()
    (app1 / ".snakebids").touch()
    (app1 / "config").mkdir()
    (app1 / "config" / "config.yaml").touch()

    (app1 / "results" / "app1").mkdir(parents=True)
    (app1 / "results" / "app1" / "result1.data").touch()
    (app1 / "results" / "app1" / "link").symlink_to(
        app1 / "workflow" / "rules" / "rule1.smk"
    )

    (app1 / "workflow" / "rules").mkdir(parents=True)
    (app1 / "workflow" / "Snakefile").touch()
    (app1 / "workflow" / "rules" / "rule1.smk").touch()

    (app1 / ".snakemake").mkdir()

    output = tmp_path / "output"
    output.mkdir()
    (output / "random_file").touch()
    (output / ".snakebids").touch()
    (output / "app1").mkdir(parents=True)
    (output / "app1" / "special_results.data").touch()

    old_style = tmp_path / "old-style"
    old_style.mkdir()
    (old_style / "config").mkdir()
    (old_style / "config" / "config.yaml").mkdir()
    (old_style / "app1").mkdir(parents=True)
    (old_style / "app1" / "special_results.data").touch()

    randomdir = tmp_path / "randomdir"
    randomdir.mkdir()
    (randomdir / "somefile").touch()

    return tmp_path


class TestCheckForResultsFolder:
    def test_raises_error_if_folder_exists(self, tmp_path: Path):
        (tmp_path / "results").mkdir()
        with pytest.raises(RunError):
            output._check_for_results_folder(tmp_path)

    def test_raises_error_if_file_exists(self, tmp_path: Path):
        (tmp_path / "results").touch()
        with pytest.raises(RunError):
            output._check_for_results_folder(tmp_path)

    def test_returns_results_path_if_doesnt_exist(self, tmp_path: Path):
        path = output._check_for_results_folder(tmp_path)
        assert path == tmp_path / "results"


class TestGetSnakebidsFile:
    def test_raises_error_if_contents_but_no_snakebids_file(self, tmp_path: Path):
        (tmp_path / "some_file").touch()
        with pytest.raises(RunError):
            output._get_snakebids_file(tmp_path)

    def test_returns_path_of_snakebids_file_if_found(self, tmp_path: Path):
        (tmp_path / "some_file").touch()
        (tmp_path / ".snakebids").touch()
        path = output._get_snakebids_file(tmp_path)
        assert path == tmp_path / ".snakebids"

    def test_returns_none_if_directory_empty(self, tmp_path: Path):
        assert output._get_snakebids_file(tmp_path) == None

    def test_returns_none_if_directory_nonexistant(self, tmp_path: Path):
        empty_dir = tmp_path / "empty"
        empty_dir.mkdir()
        assert output._get_snakebids_file(empty_dir) == None


def test_copy_snakemake_app(fake_snakemake: Path):
    # Major uncovered cases:
    #   - symlinks pointing outside the copied body
    app2 = fake_snakemake / "app2"
    out = output._copy_snakemake_app(fake_snakemake / "app1", app2)

    assert out == fake_snakemake / "app2" / "results"
    assert not (out / ".snakemake").exists()
    assert dirlen(app2 / "results") == 0
    assert dirlen(app2 / "config") == 0
    assert (app2 / "workflow/Snakefile").exists()
    assert (app2 / "workflow/rules/rule1.smk").exists()
    assert not (app2 / "results/app1/link").exists()


class TestConvertOutput:
    @pytest.mark.parametrize(
        "mode,force", it.product(("workflow", "bidsapp"), (True, False))
    )
    def test_does_nothing_when_start_and_end_are_the_same(
        self, tmp_path: Path, mode: output.Mode, force: bool
    ):
        if mode == "bidsapp":
            assert (
                output._convert_output(mode, mode, tmp_path, tmp_path, force)
                is tmp_path
            )
        else:
            assert (
                output._convert_output(mode, mode, tmp_path, tmp_path, force)
                == tmp_path / "results"
            )

    def test_workflow_conversion_fails_when_results_present(self, fake_snakemake: Path):
        with pytest.raises(RunError):
            output._convert_output(
                "bidsapp",
                "workflow",
                fake_snakemake / "app1",
                fake_snakemake / "app1",
                True,
            )

    def test_conversion_to_workflow(self, fake_snakemake: Path):
        out = output._convert_output(
            "bidsapp",
            "workflow",
            fake_snakemake / "app1",
            fake_snakemake / "output",
            True,
        )
        assert out == fake_snakemake / "output" / "results"
        assert dirlen(fake_snakemake / "output") == 5
        assert dirlen(out / "app1") == 1
        assert (out / "app1" / "special_results.data").exists()

    def test_raise_exception_when_convert_to_bidsapp_without_force(
        self, fake_snakemake: Path
    ):
        output._convert_output(
            "bidsapp",
            "workflow",
            fake_snakemake / "app1",
            fake_snakemake / "output",
            True,
        )

        with pytest.raises(RunError):
            output._convert_output(
                "workflow",
                "bidsapp",
                fake_snakemake / "app1",
                fake_snakemake / "output",
                False,
            )

    def test_converts_to_bidsapp_when_forced(self, fake_snakemake: Path):
        output._convert_output(
            "bidsapp",
            "workflow",
            fake_snakemake / "app1",
            fake_snakemake / "output",
            True,
        )

        out = output._convert_output(
            "workflow",
            "bidsapp",
            fake_snakemake / "app1",
            fake_snakemake / "output",
            True,
        )

        assert out == fake_snakemake / "output"
        assert dirlen(out) == 3


class TestPrepareOutput:
    def test_creates_new_empty_workflow_app(self, fake_snakemake: Path):
        out = output.prepare_output(
            fake_snakemake / "app1", fake_snakemake / "new-output", "workflow", False
        )

        assert out == fake_snakemake / "new-output" / "results"
        assert dirlen(fake_snakemake / "new-output") == 5
        assert out.exists()
        assert dirlen(out) == 0
        with (fake_snakemake / "new-output" / ".snakebids").open() as f:
            assert json.load(f)["mode"] == "workflow"

    def test_creates_new_empty_bidsapp(self, fake_snakemake: Path):
        out = output.prepare_output(
            fake_snakemake / "app1", fake_snakemake / "new-output", "bidsapp", False
        )

        assert out == fake_snakemake / "new-output"
        assert out.exists()
        assert dirlen(out) == 1
        with (out / ".snakebids").open() as f:
            assert json.load(f)["mode"] == "bidsapp"

    def test_converts_bidsapp_to_workflow(self, fake_snakemake: Path):
        with (fake_snakemake / "output" / ".snakebids").open("w") as f:
            json.dump({"mode": "bidsapp"}, f)
        out = output.prepare_output(
            fake_snakemake / "app1", fake_snakemake / "output", "workflow", False
        )

        assert out == fake_snakemake / "output" / "results"
        assert not (out / ".snakemake").exists()
        assert dirlen(fake_snakemake / "output") == 5
        assert out.exists()
        assert dirlen(out) == 2
        with (fake_snakemake / "output" / ".snakebids").open() as f:
            assert json.load(f)["mode"] == "workflow"

    def test_converts_workflow_app_to_bidsapp(self, fake_snakemake: Path):
        with (fake_snakemake / "output" / ".snakebids").open("w") as f:
            json.dump({"mode": "bidsapp"}, f)
        output.prepare_output(
            fake_snakemake / "app1", fake_snakemake / "output", "workflow", False
        )

        assert not (fake_snakemake / "output" / ".snakemake").exists()
        (fake_snakemake / "output" / ".snakemake").touch()

        out = output.prepare_output(
            fake_snakemake / "app1", fake_snakemake / "output", "bidsapp", True
        )

        assert out == fake_snakemake / "output"
        assert (fake_snakemake / "output" / ".snakemake").exists()
        assert out.exists()
        assert dirlen(out) == 4
        with (out / ".snakebids").open() as f:
            assert json.load(f)["mode"] == "bidsapp"

    def test_fails_when_directory_has_contents(self, fake_snakemake: Path):
        with pytest.raises(RunError):
            output.prepare_output(
                fake_snakemake / "app1", fake_snakemake / "randomdir", "workflow", False
            )


class TestRetrofitOutput:
    def test_successfully_retrofits(self, fake_snakemake: Path, mocker: MockerFixture):
        o = fake_snakemake / "old-style"
        mocker.patch("builtins.input", return_value="yes")
        out = output.retrofit_output(o, [o / "config/config.yaml"])
        assert out == True
        assert not (o / "config").exists()
        assert dirlen(o / "app1") == 1
        with (o / ".snakebids").open() as f:
            assert json.load(f)["mode"] == "bidsapp"

    def test_fails_when_user_input_not_yes(
        self, fake_snakemake: Path, mocker: MockerFixture
    ):
        o = fake_snakemake / "old-style"
        mocker.patch("builtins.input", return_value="")
        out = output.retrofit_output(o, [o / "config/config.yaml"])
        assert out == False
        assert (o / "config").exists()

    def test_fails_when_snakebids_file_exists(
        self, fake_snakemake: Path, mocker: MockerFixture
    ):
        o = fake_snakemake / "old-style"
        (o / ".snakebids").touch()
        with pytest.raises(RunError):
            output.retrofit_output(o, [o / "config/config.yaml"])

    def test_fails_if_no_config_files_provided(
        self, fake_snakemake: Path, mocker: MockerFixture
    ):
        o = fake_snakemake / "old-style"
        with pytest.raises(RunError):
            output.retrofit_output(o, [])

    def test_fails_if_unknown_files_in_config_folder(
        self, fake_snakemake: Path, mocker: MockerFixture
    ):
        o = fake_snakemake / "old-style"
        (o / "config" / "random").touch()
        with pytest.raises(RunError):
            output.retrofit_output(o, [o / "config/config.yaml"])
