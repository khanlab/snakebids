from __future__ import annotations

import json
import re
import subprocess as sp
import sys
from pathlib import Path
from typing import Any, ClassVar

import pytest

from snakebids import bidsapp
from snakebids.plugins import SnakemakeBidsApp
from snakebids.plugins import _provenance_bids as bids_prov
from snakebids.plugins import provenance as prov_mod
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


def test_explicit_instance_replaces_default_dependency():
    custom = Provenance(dataset_description=False)
    app = bidsapp.app(
        [
            SnakemakeBidsApp(Path(), configfile_path=None, snakefile_path=Path("x")),
            custom,
        ]
    )
    registered = [p for _, p in app.pm.list_name_plugin() if isinstance(p, Provenance)]
    assert len(registered) == 1
    assert registered[0] is custom


# --------------------------------------------------------------------------------------
# Edge cases and error handling
# --------------------------------------------------------------------------------------


class _FakeDist:
    def __init__(self, text: str | None):
        self.text = text

    def read_text(self, name: str) -> str | None:
        assert name == "direct_url.json"
        return self.text


def _patch_dist(monkeypatch: pytest.MonkeyPatch, direct_url: Any):
    text = direct_url if isinstance(direct_url, str | None) else json.dumps(direct_url)

    def distribution(_: str) -> _FakeDist:
        return _FakeDist(text)

    monkeypatch.setattr(bids_prov.impm, "distribution", distribution)


class TestCodeSource:
    def test_git_unavailable(self, app_repo: Path, monkeypatch: pytest.MonkeyPatch):
        def no_git(*_: Any, **__: Any):
            msg = "git"
            raise FileNotFoundError(msg)

        monkeypatch.setattr(bids_prov.sp, "run", no_git)
        assert bids_prov.git_info(app_repo) is None

    def test_repo_without_commits(self, tmp_path: Path):
        (tmp_path / "Snakefile").write_text("")
        _git(tmp_path, "init", "-q")
        _git(tmp_path, "add", "-A")
        assert bids_prov.git_info(tmp_path) is None

    def test_missing_dir(self, tmp_path: Path):
        assert bids_prov.git_info(tmp_path / "nope") is None

    def test_clean_repo_has_no_diff(self, app_repo: Path, tmp_path: Path):
        assert not bids_prov.write_git_diff(app_repo, tmp_path / "app.diff")
        assert not (tmp_path / "app.diff").exists()

    def test_large_diff_not_written(
        self, app_repo: Path, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ):
        (app_repo / "Snakefile").write_text("changed\n")
        monkeypatch.setattr(bids_prov, "_MAX_DIFF_BYTES", 10)
        assert not bids_prov.write_git_diff(app_repo, tmp_path / "app.diff")
        assert not (tmp_path / "app.diff").exists()

    def test_prefers_app_dir_checkout(self, app_repo: Path):
        source = bids_prov.code_source(app_repo, None)
        assert source is not None
        assert source["type"] == "git"

    def test_no_distribution(self, tmp_path: Path):
        assert bids_prov.code_source(tmp_path, None) is None

    def test_uninstalled_distribution(self, tmp_path: Path):
        assert bids_prov.code_source(tmp_path, "not-a-real-dist-xyz") is None

    @pytest.mark.parametrize("text", [None, "", "not json"])
    def test_unreadable_direct_url(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch, text: str | None
    ):
        _patch_dist(monkeypatch, text)
        assert bids_prov.code_source(tmp_path, "myapp") is None

    def test_editable_git_install(
        self, app_repo: Path, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ):
        _patch_dist(
            monkeypatch,
            {"url": app_repo.as_uri(), "dir_info": {"editable": True}},
        )
        source = bids_prov.code_source(tmp_path, "myapp")
        assert source is not None
        assert source["type"] == "git"
        assert source["editable"]

    def test_editable_non_git_install(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ):
        src = tmp_path / "src"
        src.mkdir()
        _patch_dist(monkeypatch, {"url": src.as_uri(), "dir_info": {"editable": True}})
        source = bids_prov.code_source(tmp_path / "other", "myapp")
        assert source == {"type": "editable", "path": str(src)}
        assert bids_prov.source_summary(source) == f"editable install at {src}"
        assert bids_prov.code_url(source) is None

    def test_vcs_install(self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch):
        _patch_dist(
            monkeypatch,
            {
                "url": "https://tok@github.com/lab/myapp.git",
                "vcs_info": {
                    "vcs": "git",
                    "commit_id": "a" * 40,
                    "requested_revision": "v1.0",
                },
            },
        )
        source = bids_prov.code_source(tmp_path, "myapp")
        assert source == {
            "type": "vcs",
            "remote": "https://github.com/lab/myapp",
            "commit": "a" * 40,
            "requested_revision": "v1.0",
        }
        assert (
            bids_prov.code_url(source)
            == f"https://github.com/lab/myapp/tree/{'a' * 40}"
        )
        assert bids_prov.source_summary(source) == (
            f"git {'a' * 12} https://github.com/lab/myapp"
        )

    def test_url_install(self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch):
        _patch_dist(monkeypatch, {"url": "https://example.org/myapp-1.0.tar.gz"})
        source = bids_prov.code_source(tmp_path, "myapp")
        url = "https://example.org/myapp-1.0.tar.gz"
        assert source == {"type": "url", "url": url}
        assert bids_prov.source_summary(source) == url

    def test_non_forge_remote_url_not_pinned(self):
        source = {"type": "git", "commit": "abc", "remote": "https://example.org/r"}
        assert bids_prov.code_url(source) == "https://example.org/r"

    def test_guess_distribution(self, tmp_path: Path):
        assert bids_prov.guess_distribution(None) is None
        assert bids_prov.guess_distribution(tmp_path / "pytest") == "pytest"
        assert bids_prov.guess_distribution(tmp_path / "not_a_package_xyz") is None


def _fake_dockerenv(monkeypatch: pytest.MonkeyPatch, *, exists: bool):
    real_exists = Path.exists

    def fake_exists(self: Path, **kwargs: Any) -> bool:
        if str(self) == "/.dockerenv":
            return exists
        return real_exists(self, **kwargs)

    monkeypatch.setattr(Path, "exists", fake_exists)


class TestDetectContainer:
    @pytest.fixture(autouse=True)
    def _clean_env(self, monkeypatch: pytest.MonkeyPatch):
        for var in [
            "APPTAINER_CONTAINER",
            "APPTAINER_NAME",
            "SINGULARITY_CONTAINER",
            "SINGULARITY_NAME",
        ]:
            monkeypatch.delenv(var, raising=False)

    def test_apptainer(self, monkeypatch: pytest.MonkeyPatch):
        monkeypatch.setenv("APPTAINER_CONTAINER", "/images/myapp_1.0.sif")
        monkeypatch.setenv("APPTAINER_NAME", "myapp_1.0.sif")
        assert bids_prov.detect_container() == {
            "Type": "apptainer",
            "Tag": "myapp_1.0.sif",
        }

    def test_singularity_uri(self, monkeypatch: pytest.MonkeyPatch):
        monkeypatch.setenv("SINGULARITY_CONTAINER", "docker://lab/myapp:1.0")
        assert bids_prov.detect_container() == {
            "Type": "singularity",
            "Tag": "myapp:1.0",
            "URI": "docker://lab/myapp:1.0",
        }

    def test_docker(self, monkeypatch: pytest.MonkeyPatch):
        _fake_dockerenv(monkeypatch, exists=True)
        assert bids_prov.detect_container() == {"Type": "docker"}

    def test_none(self, monkeypatch: pytest.MonkeyPatch):
        _fake_dockerenv(monkeypatch, exists=False)
        assert bids_prov.detect_container() is None


class TestJsonHelpers:
    def test_read_missing(self, tmp_path: Path):
        assert bids_prov.read_json(tmp_path / "x.json") is None

    @pytest.mark.parametrize("text", ["not json", "[1, 2]"])
    def test_read_invalid(self, tmp_path: Path, text: str):
        (tmp_path / "x.json").write_text(text)
        assert bids_prov.read_json(tmp_path / "x.json") is None

    def test_json_objects(self):
        assert bids_prov.json_objects("x") == []
        assert bids_prov.json_objects([{"a": 1}, "b", 2]) == [{"a": 1}]


class TestDatasetDescriptionEdgeCases:
    def test_mixed_generated_by(self, tmp_path: Path):
        (tmp_path / "dataset_description.json").write_text(
            json.dumps({"Name": "x", "GeneratedBy": ["bids::prov#a", {"Name": "b"}]})
        )
        bids_prov.update_dataset_description(tmp_path, [_rec("r1")], "myapp", None)
        gb = json.loads((tmp_path / "dataset_description.json").read_text())[
            "GeneratedBy"
        ]
        assert gb[:2] == ["bids::prov#a", {"Name": "b"}]
        assert gb[2]["Name"] == "myapp"

    def test_container_and_unknown_version(self, tmp_path: Path):
        rec = _rec("r1")
        rec.app_version = "unknown"
        rec.container = {"Type": "apptainer", "Tag": "myapp.sif"}
        (entry,) = bids_prov.generated_by_summary(
            [rec], "myapp", "code/snakebids", {"Name": "myapp", "Version": "0.0.0"}
        )
        assert "Version" not in entry
        assert entry["Container"] == rec.container

    def test_dataset_link_alias_collision(self, tmp_path: Path):
        rec = _rec("r1")
        rec.inputs = {"bids_dir": str(tmp_path / "new")}
        links = bids_prov.dataset_links([rec], {"bids_dir": "file:///other"})
        assert links == {
            "bids_dir": "file:///other",
            "bids_dir2": (tmp_path / "new").as_uri(),
        }

    def test_invalid_dataset_links_left_alone(self, tmp_path: Path):
        (tmp_path / "dataset_description.json").write_text(
            json.dumps({"Name": "x", "DatasetLinks": "bad"})
        )
        bids_prov.update_dataset_description(tmp_path, [_rec("r1")], "myapp", None)
        desc = json.loads((tmp_path / "dataset_description.json").read_text())
        assert desc["DatasetLinks"] == "bad"


class TestRecordEdgeCases:
    def test_write_requires_path(self):
        with pytest.raises(ValueError, match="no path"):
            _rec("x").write()
        with pytest.raises(ValueError, match="no path"):
            _rec("x").write_script()

    def test_script_mentions_container(self, tmp_path: Path):
        rec = _rec("x")
        rec.path = tmp_path
        rec.container = {"Type": "apptainer", "Tag": "myapp.sif"}
        rec.write_script()
        assert "apptainer container myapp.sif" in (tmp_path / "command.sh").read_text()

    def test_unreadable_record_skipped(self, tmp_path: Path):
        good = _rec("20260101T000000-aaaa")
        good.path = tmp_path / "runs" / good.run_id
        good.write()
        bad = tmp_path / "runs" / "20260102T000000-bbbb"
        bad.mkdir()
        (bad / "provenance.json").write_text("{not json")
        assert [r.run_id for r in load_history(tmp_path)] == [good.run_id]

    def test_load_history_missing_dir(self, tmp_path: Path):
        assert load_history(tmp_path) == []

    def test_prog_for_non_executable_script(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ):
        script = tmp_path / "run.py"
        script.write_text("")
        script.chmod(0o644)
        monkeypatch.setattr(sys, "argv", [str(script)])
        assert prov_mod._get_prog() == [sys.executable, str(script)]  # pyright: ignore[reportPrivateUsage]

    def test_prog_fallback(self, monkeypatch: pytest.MonkeyPatch):
        monkeypatch.setattr(sys, "argv", [""])
        assert prov_mod._get_prog() == ["app"]  # pyright: ignore[reportPrivateUsage]

    def test_safe_version_missing(self):
        assert prov_mod._safe_version("not-a-real-dist-xyz") is None  # pyright: ignore[reportPrivateUsage]


class TestPopFlag:
    def test_equals_form(self):
        assert prov_mod._pop_flag(["a", "--x=1", "b"], "--x") == ("1", ["a", "b"])  # pyright: ignore[reportPrivateUsage]

    def test_missing_value(self):
        assert prov_mod._pop_flag(["a", "--x"], "--x") == (None, ["a"])  # pyright: ignore[reportPrivateUsage]

    def test_absent(self):
        assert prov_mod._pop_flag(["a"], "--x") is None  # pyright: ignore[reportPrivateUsage]


class TestCliErrors:
    @pytest.mark.parametrize("flag", ["--history", "--replay"])
    def test_flag_requires_dir(self, flag: str, capsys: pytest.CaptureFixture[str]):
        with pytest.raises(SystemExit) as err:
            _run([flag])
        assert err.value.code == 1
        assert "requires an OUTPUT_DIR" in capsys.readouterr().err

    def test_history_empty(self, tmp_path: Path, capsys: pytest.CaptureFixture[str]):
        with pytest.raises(SystemExit) as err:
            _run(["--history", str(tmp_path)])
        assert err.value.code == 0
        assert "No recorded runs" in capsys.readouterr().out

    def test_history_shows_version_and_code(
        self,
        app_repo: Path,
        tmp_path: Path,
        monkeypatch: pytest.MonkeyPatch,
        capsys: pytest.CaptureFixture[str],
    ):
        monkeypatch.chdir(tmp_path)
        _run_bids(["out"], app_repo, tmp_path, version="2.3.4")
        _run(["--replay", "out"])
        capsys.readouterr()
        with pytest.raises(SystemExit):
            _run(["--history", "out"])
        out = capsys.readouterr().out
        assert "myapp 2.3.4" in out
        assert "code: git " in out
        assert "(replay of " in out

    def test_replay_with_missing_cwd(
        self,
        tmp_path: Path,
        monkeypatch: pytest.MonkeyPatch,
        capsys: pytest.CaptureFixture[str],
    ):
        gone = tmp_path / "gone"
        gone.mkdir()
        monkeypatch.chdir(gone)
        _run([str(tmp_path / "out")])
        monkeypatch.chdir(tmp_path)
        gone.rmdir()
        with pytest.raises(SystemExit) as err:
            _run(["--replay", str(tmp_path / "out")])
        assert err.value.code == 1
        assert "no longer exists" in capsys.readouterr().err


class _Raiser(_Recorder):
    def __init__(self, exc: BaseException):
        super().__init__()
        self.exc = exc

    @bidsapp.hookimpl
    def run(self, config: dict[str, Any]):
        raise self.exc


class TestRunStatus:
    @pytest.mark.parametrize(
        ("exc", "status", "code"),
        [
            (KeyboardInterrupt(), "interrupted", 130),
            (RuntimeError("boom"), "failed", 1),
            (SystemExit("error message"), "failed", 1),
            (SystemExit(None), "success", 0),
        ],
    )
    def test_status_recorded(
        self,
        tmp_path: Path,
        monkeypatch: pytest.MonkeyPatch,
        exc: BaseException,
        status: str,
        code: int,
    ):
        monkeypatch.chdir(tmp_path)
        with pytest.raises(type(exc)):
            bidsapp.app([Provenance(), _Raiser(exc)]).run(["out"])
        (record,) = load_history(_prov_dir(Path("out")))
        assert (record.status, record.exit_code) == (status, code)


class TestWriteFailures:
    def test_record_write_failure_is_not_fatal(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ):
        monkeypatch.chdir(tmp_path)

        def fail(*_: Any, **__: Any):
            msg = "read-only"
            raise OSError(msg)

        monkeypatch.setattr(RunRecord, "write", fail)
        rec = _run(["out"])
        assert rec.configs
        assert "plugins.provenance.run_id" not in rec.configs[0]

    def test_update_failure_is_not_fatal(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ):
        monkeypatch.chdir(tmp_path)
        real_write = RunRecord.write
        calls: list[int] = []

        def fail_second(self: RunRecord):
            calls.append(1)
            if len(calls) > 1:
                msg = "disk full"
                raise OSError(msg)
            real_write(self)

        monkeypatch.setattr(RunRecord, "write", fail_second)
        _run(["out"])
        (record,) = load_history(_prov_dir(Path("out")))
        assert record.status == "running"

    def test_dataset_description_failure_is_not_fatal(
        self, app_repo: Path, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ):
        monkeypatch.chdir(tmp_path)

        def fail(*_: Any, **__: Any):
            msg = "read-only"
            raise OSError(msg)

        monkeypatch.setattr(bids_prov, "update_dataset_description", fail)
        _run_bids(["out"], app_repo, tmp_path)
        assert load_history(_prov_dir(Path("out")))

    def test_no_user(self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch):
        monkeypatch.chdir(tmp_path)

        def fail():
            msg = "no user"
            raise KeyError(msg)

        monkeypatch.setattr(prov_mod.getpass, "getuser", fail)
        _run(["out"])
        (record,) = load_history(_prov_dir(Path("out")))
        assert record.user is None

    def test_bids_views_without_record(self):
        Provenance()._write_bids_views()  # pyright: ignore[reportPrivateUsage]

    def test_finish_without_record(self):
        Provenance()._finish("success", 0)  # pyright: ignore[reportPrivateUsage]


@pytest.mark.parametrize(
    ("seconds", "expected"),
    [
        (0, "0s"),
        (9.4, "9s"),
        (252, "4m12s"),
        (7380, "2h03m"),
        (100800, "1d04h"),
    ],
)
def test_format_duration(seconds: float, expected: str):
    assert prov_mod.format_duration(seconds) == expected


class TestDuration:
    def test_duration_from_times(self):
        rec = _rec("x")
        rec.start_time = "2026-11-06T13:40:23.5-05:00"
        rec.end_time = "2026-11-06T13:44:35.5-05:00"
        assert rec.duration == 252  # noqa: PLR2004
        assert rec.to_json()["duration_seconds"] == 252  # noqa: PLR2004

    def test_no_duration_while_running(self):
        rec = _rec("x")
        rec.start_time = "2026-11-06T13:40:23-05:00"
        assert rec.duration is None
        assert rec.to_json()["duration_seconds"] is None

    def test_invalid_times(self):
        rec = _rec("x")
        rec.end_time = "not a time"
        assert rec.duration is None

    def test_history_shows_duration(
        self,
        tmp_path: Path,
        monkeypatch: pytest.MonkeyPatch,
        capsys: pytest.CaptureFixture[str],
    ):
        monkeypatch.chdir(tmp_path)
        _run(["out"])
        (record,) = load_history(_prov_dir(Path("out")))
        assert record.duration is not None
        capsys.readouterr()
        with pytest.raises(SystemExit):
            _run(["--history", "out"])
        header = capsys.readouterr().out.splitlines()[0]
        assert re.search(r"\d{4}-\d\d-\d\d \d\d:\d\d \(\d+s\)", header)
