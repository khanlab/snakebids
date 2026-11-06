from __future__ import annotations

import json
import subprocess as sp
import sys
from pathlib import Path
from typing import Any, ClassVar

import pytest

from snakebids import bidsapp
from snakebids.plugins import _provenance_bids as bids_prov
from snakebids.plugins.provenance import (
    Provenance,
    RunRecord,
    _is_non_executing,
    _select_run,
    load_history,
)

FAIL_CODE = 3


class _Recorder:
    """Plugin capturing the config passed to run and optionally failing."""

    def __init__(self, exit_code: int | None = None):
        self.exit_code = exit_code
        self.configs: list[dict[str, Any]] = []

    @bidsapp.hookimpl
    def add_cli_arguments(self, parser: Any):
        parser.add_argument("output_dir", type=Path)
        parser.add_argument("--opt", default="a")

    @bidsapp.hookimpl
    def run(self, config: dict[str, Any]):
        self.configs.append(dict(config))
        if self.exit_code is not None:
            sys.exit(self.exit_code)


def _run(args: list[str], exit_code: int | None = None) -> _Recorder:
    rec = _Recorder(exit_code)
    bidsapp.app([Provenance(), rec]).run(args)
    return rec


def _prov_dir(out: Path) -> Path:
    return out / "code" / "snakebids"


class TestRecording:
    def test_run_is_recorded(self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch):
        monkeypatch.chdir(tmp_path)
        rec = _run(["out", "--opt", "b"])
        (record,) = load_history(_prov_dir(Path("out")))
        assert record.argv == ["out", "--opt", "b"]
        assert record.cwd == str(tmp_path)
        assert record.status == "success"
        assert record.exit_code == 0
        assert rec.configs[0]["plugins.provenance.run_id"] == record.run_id
        assert record.path is not None
        assert (record.path / "config.yaml").exists()
        script = (record.path / "command.sh").read_text()
        assert "--opt b" in script
        assert f"cd {tmp_path}" in script

    def test_failure_status_recorded(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ):
        monkeypatch.chdir(tmp_path)
        with pytest.raises(SystemExit):
            _run(["out"], exit_code=FAIL_CODE)
        (record,) = load_history(_prov_dir(Path("out")))
        assert record.status == "failed"
        assert record.exit_code == FAIL_CODE

    def test_dry_runs_not_recorded(self, tmp_path: Path):
        plugin = Provenance()
        plugin.get_argv(["x"], {})
        plugin.finalize_config({"output_dir": tmp_path, "snakemake_args": ["-np"]})
        assert load_history(_prov_dir(tmp_path)) == []

    def test_no_output_dir_does_nothing(self):
        plugin = Provenance()
        plugin.finalize_config({})
        assert plugin._record is None  # pyright: ignore[reportPrivateUsage]


class TestReplay:
    def test_replay_reuses_args_and_appends_extra(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ):
        monkeypatch.chdir(tmp_path)
        _run(["out", "--opt", "b"])
        _run(["out", "--opt", "c"])
        elsewhere = tmp_path / "elsewhere"
        elsewhere.mkdir()
        monkeypatch.chdir(elsewhere)

        rec = _run(["--replay", str(tmp_path / "out") + "@1", "--opt", "z"])
        assert rec.configs[0]["opt"] == "z"
        # chdir'd back to original cwd so relative output_dir resolves the same
        assert Path.cwd() == tmp_path

        history = load_history(_prov_dir(tmp_path / "out"))
        assert len(history) == 3  # noqa: PLR2004
        assert history[-1].replayed_from == history[0].run_id
        assert history[-1].argv == ["out", "--opt", "b", "--opt", "z"]

    def test_replay_missing_history_exits(
        self, tmp_path: Path, capsys: pytest.CaptureFixture[str]
    ):
        with pytest.raises(SystemExit) as err:
            _run(["--replay", str(tmp_path)])
        assert err.value.code == 1
        assert "No recorded runs" in capsys.readouterr().err

    def test_history_prints_and_exits(
        self,
        tmp_path: Path,
        monkeypatch: pytest.MonkeyPatch,
        capsys: pytest.CaptureFixture[str],
    ):
        monkeypatch.chdir(tmp_path)
        _run(["out", "--opt", "b"])
        capsys.readouterr()
        with pytest.raises(SystemExit) as err:
            _run(["--history", "out"])
        assert err.value.code == 0
        out = capsys.readouterr().out
        assert "[1]" in out
        assert "--opt b" in out


def _rec(run_id: str) -> RunRecord:
    return RunRecord(run_id=run_id, argv=[], cwd=".", prog=["app"], start_time="")


class TestSelectRun:
    records: ClassVar[list[RunRecord]] = [
        _rec("20260101T000000-aaaa"),
        _rec("20260102T000000-bbbb"),
    ]

    @pytest.mark.parametrize(
        ("selector", "expected"),
        [(None, 1), ("latest", 1), ("1", 0), ("-1", 1), ("-2", 0), ("20260101", 0)],
    )
    def test_selectors(self, selector: str | None, expected: int):
        assert _select_run(self.records, selector) is self.records[expected]

    @pytest.mark.parametrize("selector", ["0", "3", "-3", "2026", "nope"])
    def test_bad_selectors(self, selector: str):
        with pytest.raises(LookupError):
            _select_run(self.records, selector)


@pytest.mark.parametrize(
    ("args", "expected"),
    [
        ([], False),
        (["-c4"], False),
        (["--cores", "4", "-p"], False),
        (["-n"], True),
        (["-np"], True),
        (["--dry-run"], True),
        (["--unlock"], True),
    ],
)
def test_is_non_executing(args: list[str], expected: bool):
    assert _is_non_executing(args) is expected


def test_record_round_trip(tmp_path: Path):
    rec = _rec("x")
    rec.path = tmp_path / "x"
    rec.write()
    data = json.loads((tmp_path / "x" / "provenance.json").read_text())
    assert data["command"] == "app"
    assert RunRecord.from_json(data) == rec


# --------------------------------------------------------------------------------------
# BIDS views: dataset_description.json, prov/, git source
# --------------------------------------------------------------------------------------


def _git(repo: Path, *args: str) -> None:
    sp.run(
        ["git", "-c", "user.name=t", "-c", "user.email=t@t", *args],
        cwd=repo,
        check=True,
        capture_output=True,
    )


@pytest.fixture
def app_repo(tmp_path: Path) -> Path:
    repo = tmp_path / "myapp"
    repo.mkdir()
    (repo / "Snakefile").write_text("rule all:\n")
    (repo / "pipeline_description.json").write_text(
        json.dumps(
            {
                "Name": "Dataset generated by myapp",
                "BIDSVersion": "1.8.0",
                "DatasetType": "derivative",
                "GeneratedBy": [
                    {"Name": "myapp", "Version": "0.0.0", "Author": "Someone"}
                ],
            }
        )
    )
    _git(repo, "init", "-q")
    _git(repo, "add", "-A")
    _git(repo, "commit", "-qm", "init")
    _git(repo, "remote", "add", "origin", "https://me:SECRET@github.com/lab/myapp.git")
    return repo


class _BidsRecorder(_Recorder):
    def __init__(self, app_dir: Path, bids_dir: Path, version: str | None):
        super().__init__()
        self.app_dir = app_dir
        self.bids_dir = bids_dir
        self.version = version

    @bidsapp.hookimpl
    def initialize_config(self, config: dict[str, Any]):
        config["snakemake_dir"] = self.app_dir
        config["bids_dir"] = self.bids_dir
        if self.version:
            config["plugins.version.version"] = self.version


def _run_bids(
    args: list[str],
    app_dir: Path,
    bids_dir: Path,
    version: str | None = "1.0.0",
    **kwargs: Any,
):
    recorder = _BidsRecorder(app_dir, bids_dir, version)
    bidsapp.app([Provenance(**kwargs), recorder]).run(args)


@pytest.mark.parametrize(
    ("url", "expected"),
    [
        ("git@github.com:lab/app.git", "https://github.com/lab/app"),
        ("https://user:tok@github.com/lab/app.git", "https://github.com/lab/app"),
        ("ssh://git@host.org/lab/app.git", "https://host.org/lab/app"),
        ("https://github.com/lab/app", "https://github.com/lab/app"),
    ],
)
def test_sanitize_url(url: str, expected: str):
    assert bids_prov.sanitize_url(url) == expected


class TestGitSource:
    def test_clean_checkout(self, app_repo: Path):
        info = bids_prov.git_info(app_repo)
        assert info is not None
        assert not info["dirty"]
        assert info["remote"] == "https://github.com/lab/myapp"
        assert len(info["commit"]) == 40  # noqa: PLR2004

    def test_untracked_dir_is_ignored(self, app_repo: Path):
        venv = app_repo / ".venv" / "myapp"
        venv.mkdir(parents=True)
        (app_repo / ".gitignore").write_text(".venv\n")
        assert bids_prov.git_info(venv) is None

    def test_not_a_repo(self, tmp_path: Path):
        assert bids_prov.git_info(tmp_path) is None

    def test_dirty_checkout_saves_diff(
        self, app_repo: Path, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ):
        monkeypatch.chdir(tmp_path)
        (app_repo / "Snakefile").write_text("rule all:\n  # changed\n")
        _run_bids(["out"], app_repo, tmp_path)
        (record,) = load_history(_prov_dir(Path("out")))
        assert record.source is not None
        assert record.source["dirty"]
        assert record.path is not None
        assert "# changed" in (record.path / "app.diff").read_text()


class TestDatasetDescription:
    def test_created_from_template(
        self, app_repo: Path, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ):
        monkeypatch.chdir(tmp_path)
        _run_bids(["out", "--opt", "x"], app_repo, tmp_path)
        _run_bids(["out", "--opt", "y"], app_repo, tmp_path)
        text = Path("out/dataset_description.json").read_text()
        desc = json.loads(text)
        assert desc["Name"] == "Dataset generated by myapp"
        assert desc["BIDSVersion"] == "1.8.0"
        # template entry replaced by a single summary entry inheriting its keys
        (entry,) = desc["GeneratedBy"]
        assert entry["Name"] == "myapp"
        assert entry["Author"] == "Someone"
        assert entry["Version"] == "1.0.0"
        assert entry["Description"].startswith("Recorded by snakebids: 2 runs")
        assert "code/snakebids/" in entry["Description"]
        assert entry["CodeURL"].startswith("https://github.com/lab/myapp/tree/")
        # summary only: no per-run commands
        assert "--opt" not in text
        assert "SECRET" not in text
        assert desc["DatasetLinks"] == {"bids_dir": tmp_path.resolve().as_uri()}

    def test_one_entry_per_version(
        self, app_repo: Path, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ):
        monkeypatch.chdir(tmp_path)
        for version in ["1.0.0", "1.0.0", "1.1.0", None]:
            _run_bids(["out"], app_repo, tmp_path, version=version)
        gb = json.loads(Path("out/dataset_description.json").read_text())["GeneratedBy"]
        assert [e.get("Version") for e in gb] == ["1.0.0", "1.1.0", None]
        assert gb[0]["Description"].startswith("Recorded by snakebids: 2 runs")
        assert gb[1]["Description"].startswith("Recorded by snakebids: 1 run,")

    def test_failed_runs_counted(
        self, app_repo: Path, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ):
        monkeypatch.chdir(tmp_path)
        _run_bids(["out"], app_repo, tmp_path)
        out = Path("out")
        recorder = _BidsRecorder(app_repo, tmp_path, "1.0.0")
        recorder.exit_code = 1
        with pytest.raises(SystemExit):
            bidsapp.app([Provenance(), recorder]).run(["out"])
        (entry,) = json.loads((out / "dataset_description.json").read_text())[
            "GeneratedBy"
        ]
        assert "2 runs (1 failed)" in entry["Description"]

    def test_foreign_entries_and_keys_preserved(
        self, app_repo: Path, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ):
        monkeypatch.chdir(tmp_path)
        out = Path("out")
        out.mkdir()
        (out / "dataset_description.json").write_text(
            json.dumps(
                {
                    "Name": "custom",
                    "BIDSVersion": "1.9.0",
                    "License": "CC0",
                    "GeneratedBy": [{"Name": "othertool", "Version": "1"}],
                }
            )
        )
        _run_bids(["out"], app_repo, tmp_path)
        _run_bids(["out"], app_repo, tmp_path)
        desc = json.loads((out / "dataset_description.json").read_text())
        assert desc["Name"] == "custom"
        assert desc["License"] == "CC0"
        names = [e["Name"] for e in desc["GeneratedBy"]]
        assert names == ["othertool", "myapp"]

    def test_generated_by_id_form_left_alone(
        self, app_repo: Path, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ):
        monkeypatch.chdir(tmp_path)
        out = Path("out")
        out.mkdir()
        (out / "dataset_description.json").write_text(
            json.dumps({"Name": "x", "GeneratedBy": ["bids::prov#other-1234"]})
        )
        _run_bids(["out"], app_repo, tmp_path)
        gb = json.loads((out / "dataset_description.json").read_text())["GeneratedBy"]
        assert gb == ["bids::prov#other-1234"]

    def test_disabled(
        self, app_repo: Path, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ):
        monkeypatch.chdir(tmp_path)
        _run_bids(["out"], app_repo, tmp_path, dataset_description=False)
        assert not Path("out/dataset_description.json").exists()
        assert load_history(_prov_dir(Path("out")))

    def test_no_prov_dir(
        self, app_repo: Path, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ):
        monkeypatch.chdir(tmp_path)
        _run_bids(["out"], app_repo, tmp_path)
        assert not Path("out/prov").exists()


def test_history_omits_unknown_version(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
):
    monkeypatch.chdir(tmp_path)
    plugin = Provenance()
    plugin.get_argv(["out"], {})
    plugin.finalize_config(
        {"output_dir": Path("out"), "plugins.version.version": "unknown"}
    )
    with pytest.raises(SystemExit):
        _run(["--history", "out"])
    header = capsys.readouterr().out.splitlines()[0]
    assert "unknown" not in header
