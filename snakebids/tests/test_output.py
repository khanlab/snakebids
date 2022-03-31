from __future__ import absolute_import

import itertools as it
import json
from pathlib import Path
from typing import Callable

import pytest
from pytest_mock.plugin import MockerFixture

import snakebids.utils.output as output
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


class TestGetSnakebidsFile:
    def test_raises_error_if_contents_but_no_snakebids_file(self, tmp_path: Path):
        (tmp_path / "some_file").touch()
        with pytest.raises(RunError):
            output._get_snakebids_file(tmp_path)

    def test_raises_error_if_malformed_snakebids_file_found(self, tmp_path: Path):
        (tmp_path / "some_file").touch()
        (tmp_path / ".snakebids").touch()
        with pytest.raises(RunError):
            output._get_snakebids_file(tmp_path)

    def test_returns_contents_of_snakebids_file_if_valid(self, tmp_path: Path):
        (tmp_path / "some_file").touch()
        (tmp_path / ".snakebids").touch()
        with (tmp_path / ".snakebids").open("w") as f:
            json.dump({"mode": "foo"}, f)
        path = output._get_snakebids_file(tmp_path)
        assert path == {"mode": "foo"}

    def test_returns_none_if_directory_empty(self, tmp_path: Path):
        assert output._get_snakebids_file(tmp_path) is None

    def test_returns_none_if_directory_nonexistant(self, tmp_path: Path):
        empty_dir = tmp_path / "empty"
        empty_dir.mkdir()
        assert output._get_snakebids_file(empty_dir) is None


class TestPrepareBidsappOutput:
    def test_creates_new_empty_bidsapp(self, fake_snakemake: Path):
        out = fake_snakemake / "new-output"
        output.prepare_bidsapp_output(out, False)

        assert out.exists()
        assert dirlen(out) == 1
        with (out / ".snakebids").open() as f:
            assert json.load(f)["mode"] == "bidsapp"

    def test_fails_when_directory_has_contents(self, fake_snakemake: Path):
        with pytest.raises(RunError):
            output.prepare_bidsapp_output(fake_snakemake / "randomdir", False)

    def test_success_on_directory_with_contents_when_forced(self, fake_snakemake: Path):
        out = fake_snakemake / "randomdir"
        output.prepare_bidsapp_output(out, True)
        assert out.exists()
        assert dirlen(out) == 2
        with (out / ".snakebids").open() as f:
            assert json.load(f)["mode"] == "bidsapp"
