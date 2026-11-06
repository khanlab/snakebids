"""Record the CLI invocations that produced an output dataset, and replay them.

Every (non-dry) run of the app writes a run record into the output directory::

    <output_dir>/code/snakebids/runs/<run_id>/
        command.sh        # re-runnable shell script (cd + original command)
        config.yaml       # snapshot of the fully resolved config for this run
        provenance.json   # argv, cwd, versions, git commit, host, container, status
        app.diff          # uncommitted changes, if the app is a dirty git checkout

``code/`` is the BIDS-sanctioned location for code used to generate a dataset, so it
is ignored by the BIDS validator and travels with the dataset.

The records can then be inspected or replayed from the CLI::

    myapp --history OUTPUT_DIR              # list runs and their commands
    myapp --replay OUTPUT_DIR               # re-run the latest run
    myapp --replay OUTPUT_DIR@2 --cores 8   # re-run run #2 with extra/overriding args
"""

from __future__ import annotations

import argparse
import getpass
import importlib.metadata as impm
import json
import logging
import os
import platform
import re
import secrets
import shlex
import sys
from collections.abc import Generator, Sequence
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, NoReturn

import attrs

from snakebids import bidsapp
from snakebids.io.config import write_config
from snakebids.plugins import _provenance_bids as bids_prov
from snakebids.plugins.base import PluginBase

__all__ = ["Provenance", "RunRecord", "load_history"]

logger = logging.getLogger(__name__)

_RECORD_FILE = "provenance.json"

#: Snakemake flags that mean no outputs will actually be produced
_NON_EXECUTING_FLAGS = frozenset(
    {
        "-n",
        "--dry-run",
        "--dryrun",
        "--dry_run",
        "--dag",
        "--rulegraph",
        "--filegraph",
        "--d3dag",
        "--lint",
        "-l",
        "--list",
        "--list-rules",
        "--list-target-rules",
        "--lt",
        "-S",
        "--summary",
        "-D",
        "--detailed-summary",
        "--unlock",
        "--cleanup-metadata",
        "--list-changes",
        "--lc",
        "--list-input-changes",
        "--li",
        "--list-params-changes",
        "--lp",
        "--list-code-changes",
        "--list-untracked",
        "--lu",
        "--archive",
        "--containerize",
        "--generate-unit-tests",
        "--help",
        "-h",
    }
)


def _is_non_executing(snakemake_args: Sequence[str]) -> bool:
    for arg in snakemake_args:
        if arg in _NON_EXECUTING_FLAGS:
            return True
        # clustered short flags, e.g. ``-np``
        if re.fullmatch(r"-[a-zA-Z]{2,}", arg) and "n" in arg[1:]:
            return True
    return False


def _now() -> datetime:
    return datetime.now(timezone.utc).astimezone()


def _safe_version(dist: str) -> str | None:
    try:
        return impm.version(dist)
    except impm.PackageNotFoundError:
        return None


def format_duration(seconds: float) -> str:
    """Format a duration compactly, e.g. ``9s``, ``4m12s``, ``2h03m``, ``1d04h``."""
    total = round(seconds)
    days, rem = divmod(total, 86400)
    hours, rem = divmod(rem, 3600)
    minutes, secs = divmod(rem, 60)
    if days:
        return f"{days}d{hours:02d}h"
    if hours:
        return f"{hours}h{minutes:02d}m"
    if minutes:
        return f"{minutes}m{secs:02d}s"
    return f"{secs}s"


def _get_prog() -> list[str]:
    """Return the command used to launch the app, made re-runnable if possible."""
    argv0 = sys.argv[0] if sys.argv and sys.argv[0] else "app"
    path = Path(argv0)
    # ``python run.py`` or ``python -m app``: the script itself may not be executable
    if path.suffix == ".py" and path.is_file() and not os.access(path, os.X_OK):
        return [sys.executable, argv0]
    return [argv0]


@attrs.define
class RunRecord:
    """A single recorded invocation of the app."""

    run_id: str
    argv: list[str]
    cwd: str
    prog: list[str]
    start_time: str
    status: str = "running"
    end_time: str | None = None
    exit_code: int | None = None
    replayed_from: str | None = None
    app_version: str | None = None
    snakebids_version: str | None = None
    snakemake_version: str | None = None
    app_name: str | None = None
    inputs: dict[str, str] = attrs.field(factory=dict[str, str])
    """Input datasets, keyed by config key (e.g. ``bids_dir``)"""
    source: dict[str, Any] | None = None
    """Where the app code came from (git commit, pip VCS install, ...)"""
    container: dict[str, Any] | None = None
    python: str = attrs.field(factory=lambda: sys.executable)
    python_version: str = attrs.field(factory=platform.python_version)
    operating_system: str = attrs.field(factory=platform.platform)
    host: str = attrs.field(factory=platform.node)
    user: str | None = None
    path: Path | None = attrs.field(default=None, eq=False)

    @property
    def duration(self) -> float | None:
        """Run time in seconds, or None if the run hasn't finished."""
        if self.end_time is None:
            return None
        try:
            start = datetime.fromisoformat(self.start_time)
            end = datetime.fromisoformat(self.end_time)
        except ValueError:
            return None
        return (end - start).total_seconds()

    @property
    def command(self) -> str:
        """Shell-quoted command line, as it would be typed."""
        return shlex.join([*self.prog, *self.argv])

    def to_json(self) -> dict[str, Any]:
        """Serialize to a JSON-compatible dict."""
        data = attrs.asdict(self, filter=lambda a, _: a.name != "path")
        data["command"] = self.command
        data["duration_seconds"] = (
            round(self.duration, 3) if self.duration is not None else None
        )
        return data

    @classmethod
    def from_json(cls, data: dict[str, Any], path: Path | None = None) -> RunRecord:
        """Load from a dict produced by :meth:`to_json`, ignoring unknown keys."""
        fields = {a.name for a in attrs.fields(cls)} - {"path"}
        return cls(**{k: v for k, v in data.items() if k in fields}, path=path)

    def write(self) -> None:
        """Atomically write ``provenance.json`` into :attr:`path`."""
        if self.path is None:
            msg = "RunRecord has no path"
            raise ValueError(msg)
        self.path.mkdir(parents=True, exist_ok=True)
        tmp = self.path / f".{_RECORD_FILE}.tmp"
        tmp.write_text(json.dumps(self.to_json(), indent=2) + "\n")
        tmp.replace(self.path / _RECORD_FILE)

    def write_script(self) -> None:
        """Write an executable ``command.sh`` reproducing this run."""
        if self.path is None:
            msg = "RunRecord has no path"
            raise ValueError(msg)
        lines = [
            "#!/usr/bin/env bash",
            f"# Recorded by snakebids {self.snakebids_version} (run {self.run_id})",
            f"# started {self.start_time[:19]} on {self.host}"
            + (f" by {self.user}" if self.user else ""),
            f"# app version: {self.app_version}",
        ]
        if summary := bids_prov.source_summary(self.source):
            lines.append(f"# app code: {summary}")
        if self.container:
            kind, tag = self.container.get("Type"), self.container.get("Tag", "")
            lines.append(f"# NOTE: originally run inside {kind} container {tag}")
        lines += [
            "set -euo pipefail",
            f"cd {shlex.quote(self.cwd)}",
            self.command + ' "$@"',
            "",
        ]
        script = self.path / "command.sh"
        script.write_text("\n".join(lines))
        script.chmod(0o755)


def load_history(provenance_dir: Path) -> list[RunRecord]:
    """Load all run records under ``provenance_dir``, oldest first."""
    runs_dir = provenance_dir / "runs"
    if not runs_dir.is_dir():
        return []
    records: list[RunRecord] = []
    for record_file in sorted(runs_dir.glob(f"*/{_RECORD_FILE}")):
        try:
            data = json.loads(record_file.read_text())
        except (OSError, json.JSONDecodeError) as err:
            logger.warning("Skipping unreadable run record %s: %s", record_file, err)
            continue
        records.append(RunRecord.from_json(data, path=record_file.parent))
    # run_ids only have second resolution, so order by the precise start time
    return sorted(records, key=lambda r: (r.start_time, r.run_id))


def _select_run(records: list[RunRecord], selector: str | None) -> RunRecord:
    """Select a run by 1-based index (negative counts from the end) or run_id."""
    if not records:
        msg = "No recorded runs found"
        raise LookupError(msg)
    if selector is None or selector == "latest":
        return records[-1]
    # short integers are indices; anything longer is treated as a run_id prefix
    if re.fullmatch(r"-?\d{1,4}", selector):
        idx = int(selector)
        if idx == 0 or abs(idx) > len(records):
            msg = f"Run index {idx} out of range (1-{len(records)})"
            raise LookupError(msg)
        return records[idx - 1 if idx > 0 else idx]
    matches = [r for r in records if r.run_id.startswith(selector)]
    if len(matches) == 1:
        return matches[0]
    msg = (
        f"Run '{selector}' is ambiguous" if matches else f"No run matching '{selector}'"
    )
    raise LookupError(msg)


def _pop_flag(argv: list[str], flag: str) -> tuple[str | None, list[str]] | None:
    """Remove ``flag VALUE`` or ``flag=VALUE`` from argv, returning value and rest."""
    for i, arg in enumerate(argv):
        if arg == flag:
            if i + 1 >= len(argv):
                return None, argv[:i]
            return argv[i + 1], argv[:i] + argv[i + 2 :]
        if arg.startswith(f"{flag}="):
            return arg.split("=", 1)[1], argv[:i] + argv[i + 1 :]
    return None


@attrs.define(kw_only=True)
class Provenance(PluginBase):
    """Record the command line used for each run, with the option to replay it.

    On every run that will actually execute the workflow, a record is written to
    ``<output_dir>/<directory>/runs/<run_id>/``, containing the exact arguments,
    working directory, a snapshot of the resolved config, software versions, the git
    commit of the app (with a diff of any uncommitted changes), host, container, and
    the final exit status. A summary is added to the output's
    ``dataset_description.json``.

    This plugin is a dependency of :class:`~snakebids.plugins.SnakemakeBidsApp`, so
    it is enabled by default in snakemake apps. To configure it, pass it explicitly
    alongside ``SnakemakeBidsApp``. See :ref:`run-history` for a user-facing
    overview.

    Adds the following CLI flags (handled before normal argument parsing, so the
    usual required positional arguments may be omitted):

    ``--history OUTPUT_DIR``
        Print all recorded runs for an output directory and exit.
    ``--replay OUTPUT_DIR[@RUN]``
        Re-run a recorded command. ``RUN`` is a 1-based index (negative values count
        from the latest), a run_id prefix, or ``latest`` (default). Any additional
        arguments are appended to the recorded ones, so they can be used to add or
        override options (e.g. ``--cores 8``). The app first changes to the original
        working directory, so relative paths resolve as they did originally.

    Parameters
    ----------
    directory
        Location of the provenance records, relative to the output directory.
    record_dry_runs
        If True, also record runs that do not execute the workflow (e.g. ``-n``,
        ``--dag``, ``--unlock``).
    argument_group
        Title of the CLI argument group for the provenance flags.
    name
        App name used in ``dataset_description.json``. Defaults to the
        distribution name, then the ``Name`` of the first ``GeneratedBy`` entry in
        ``pipeline_description.json``, then the name of the app directory.
    distribution
        Name of the app's installed distribution, used to find out how it was
        installed (e.g. ``pip install git+https://...`` records the commit). Defaults
        to the distribution providing the top-level package named like ``app_dir``.
    app_dir
        Root of the app's code, checked for a git checkout. Defaults to the
        ``snakemake_dir`` set by :class:`~snakebids.plugins.SnakemakeBidsApp`.
    pipeline_description
        Path (relative to ``app_dir``) of a template for ``dataset_description.json``,
        used when the output dataset doesn't have one yet.
    dataset_description
        Create/update ``dataset_description.json`` in the output directory. Its
        ``GeneratedBy`` holds a compact summary of the runs (one entry per app
        version, with the latest commit and run count); full details stay in
        ``directory``. It is regenerated from the run records at the start and end
        of every run.
    """

    directory: str = "code/snakebids"
    record_dry_runs: bool = False
    argument_group: str | None = None
    name: str | None = None
    distribution: str | None = None
    app_dir: Path | None = attrs.field(
        default=None, converter=attrs.converters.optional(Path)
    )
    pipeline_description: str = "pipeline_description.json"
    dataset_description: bool = True

    PREFIX = "plugins.provenance"

    _argv: list[str] = attrs.field(factory=list[str], init=False)
    _cwd: str = attrs.field(default="", init=False)
    _replayed_from: str | None = attrs.field(default=None, init=False)
    _record: RunRecord | None = attrs.field(default=None, init=False)
    _output_dir: Path | None = attrs.field(default=None, init=False)
    _template: dict[str, Any] | None = attrs.field(default=None, init=False)

    def __eq__(self, other: Any):
        # Only one instance per app: lets an explicitly configured Provenance replace
        # the default one loaded as a dependency of SnakemakeBidsApp
        return isinstance(other, self.__class__)

    def _provenance_dir(self, output_dir: str | os.PathLike[str]) -> Path:
        return Path(output_dir) / self.directory

    @staticmethod
    def _split_selector(target: str) -> tuple[Path, str | None]:
        if "@" in target and not Path(target).exists():
            path, selector = target.rsplit("@", 1)
            return Path(path), selector
        return Path(target), None

    @bidsapp.hookimpl
    def add_cli_arguments(
        self, parser: argparse.ArgumentParser, argument_groups: Any
    ) -> None:
        """Add (documentation-only) entries for --history and --replay."""
        group = (
            argument_groups[self.argument_group]
            if self.argument_group
            else parser.add_argument_group(
                "PROVENANCE",
                "Each run is recorded in the output directory under "
                f"'{self.directory}'. These flags may be given in place of the "
                "usual positional arguments.",
            )
        )
        self.try_add_argument(
            group,
            "--history",
            dest="history",
            metavar="OUTPUT_DIR",
            help="List the recorded runs (and their commands) for OUTPUT_DIR and exit",
        )
        self.try_add_argument(
            group,
            "--replay",
            dest="replay",
            metavar="OUTPUT_DIR[@RUN]",
            help=(
                "Re-run a recorded command for OUTPUT_DIR. RUN is a 1-based index "
                "(negative counts back from latest), a run id, or 'latest' (default). "
                "Extra arguments are appended to the recorded command, e.g. to "
                "override --cores."
            ),
        )

    @bidsapp.hookimpl(tryfirst=True)
    def get_argv(self, argv: list[str], config: dict[str, Any]) -> list[str] | None:
        """Intercept --history and --replay before the parser sees the arguments."""
        self._argv = list(argv)
        self._cwd = os.getcwd()

        if (popped := _pop_flag(argv, "--history")) is not None:
            target, _ = popped
            if target is None:
                self._exit("--history requires an OUTPUT_DIR")
            self._print_history(Path(target))
            sys.exit(0)

        if (popped := _pop_flag(argv, "--replay")) is None:
            return None

        target, extra = popped
        if target is None:
            self._exit("--replay requires an OUTPUT_DIR")
        output_dir, selector = self._split_selector(target)
        try:
            run = _select_run(load_history(self._provenance_dir(output_dir)), selector)
        except LookupError as err:
            self._exit(f"Cannot replay from '{output_dir}': {err}")

        if Path(run.cwd).resolve() != Path.cwd().resolve():
            if not Path(run.cwd).is_dir():
                self._exit(
                    f"Original working directory '{run.cwd}' no longer exists; "
                    f"run {run.path / 'command.sh'} manually after editing paths"
                    if run.path
                    else f"Original working directory '{run.cwd}' no longer exists"
                )
            logger.info("Changing to original working directory: %s", run.cwd)
            os.chdir(run.cwd)
            self._cwd = run.cwd

        new_argv = [*run.argv, *extra]
        self._argv = new_argv
        self._replayed_from = run.run_id
        print(f"Replaying run {run.run_id}:\n  {shlex.join([*run.prog, *new_argv])}")
        return new_argv

    @bidsapp.hookimpl
    def update_cli_namespace(self, namespace: dict[str, Any]) -> None:
        """Drop the documentation-only arguments from the namespace."""
        self.pop(namespace, "history", None)
        self.pop(namespace, "replay", None)

    @bidsapp.hookimpl(trylast=True)
    def finalize_config(self, config: dict[str, Any]) -> None:
        """Write the run record, after all other plugins have finalized config."""
        if "output_dir" not in config:
            logger.debug("No output_dir in config; provenance will not be recorded")
            return
        if not self.record_dry_runs and _is_non_executing(
            config.get("snakemake_args", [])
        ):
            return

        start = _now()
        run_id = f"{start:%Y%m%dT%H%M%S}-{secrets.token_hex(2)}"
        output_dir = Path(config["output_dir"])
        prov_dir = self._provenance_dir(output_dir)
        try:
            user = getpass.getuser()
        except Exception:  # noqa: BLE001
            user = None

        app_dir = self.app_dir or (
            Path(config["snakemake_dir"]) if config.get("snakemake_dir") else None
        )
        distribution = self.distribution or bids_prov.guess_distribution(app_dir)
        self._template = self._load_template(app_dir)
        name = (
            self.name
            or distribution
            or self._template_name()
            or (app_dir.name if app_dir else Path(sys.argv[0]).stem)
        )
        inputs = {
            key: str(Path(config[key]).resolve())
            for key in ("bids_dir",)
            if config.get(key)
        }

        record = RunRecord(
            run_id=run_id,
            argv=self._argv,
            cwd=self._cwd,
            prog=_get_prog(),
            start_time=start.isoformat(),
            replayed_from=self._replayed_from,
            app_version=config.get("plugins.version.version"),
            snakebids_version=_safe_version("snakebids"),
            snakemake_version=_safe_version("snakemake"),
            user=user,
            app_name=name,
            inputs=inputs,
            source=bids_prov.code_source(app_dir, distribution),
            container=bids_prov.detect_container(),
            path=prov_dir / "runs" / run_id,
        )
        try:
            record.write()
            record.write_script()
            if record.source and record.source.get("dirty") and app_dir:
                bids_prov.write_git_diff(
                    Path(record.source.get("toplevel", app_dir)),
                    record.path / "app.diff",  # type: ignore[operator]
                )
            write_config(record.path / "config.yaml", config, force_overwrite=True)  # type: ignore[operator]
        except OSError as err:
            logger.warning("Could not write provenance record: %s", err)
            return
        self._record = record
        self._output_dir = output_dir
        self.assign(config, "run_id", run_id)
        self._write_bids_views()

    def _load_template(self, app_dir: Path | None) -> dict[str, Any] | None:
        if app_dir is None:
            return None
        path = app_dir / self.pipeline_description
        return bids_prov.read_json(path) if path.exists() else None

    def _template_name(self) -> str | None:
        for entry in bids_prov.json_objects((self._template or {}).get("GeneratedBy")):
            if isinstance(name := entry.get("Name"), str):
                return name
        return None

    def _write_bids_views(self) -> None:
        """Regenerate dataset_description.json from all run records."""
        if self._output_dir is None or self._record is None:
            return
        if not self.dataset_description:
            return
        out = self._output_dir
        name = self._record.app_name or "snakebids_app"
        # only include runs of this app, in case the directory is shared
        records = [
            r
            for r in load_history(self._provenance_dir(out))
            if (r.app_name or name) == name
        ]
        try:
            bids_prov.update_dataset_description(
                out, records, name, self._template, self.directory
            )
        except OSError as err:
            logger.warning("Could not write BIDS provenance: %s", err)

    @bidsapp.hookimpl(wrapper=True)
    def run(self, config: dict[str, Any]) -> Generator[None, None, None]:
        """Record the end time and exit status of the run."""
        try:
            yield
        except SystemExit as err:
            code = (
                err.code
                if isinstance(err.code, int)
                else (0 if err.code is None else 1)
            )
            self._finish("success" if code == 0 else "failed", code)
            raise
        except KeyboardInterrupt:
            self._finish("interrupted", 130)
            raise
        except BaseException:
            self._finish("failed", 1)
            raise
        self._finish("success", 0)

    def _finish(self, status: str, code: int) -> None:
        if self._record is None:
            return
        self._record.status = status
        self._record.exit_code = code
        self._record.end_time = _now().isoformat()
        try:
            self._record.write()
        except OSError as err:
            logger.warning("Could not update provenance record: %s", err)
        self._write_bids_views()

    def _print_history(self, output_dir: Path) -> None:
        records = load_history(self._provenance_dir(output_dir))
        if not records:
            print(f"No recorded runs in {self._provenance_dir(output_dir)}")
            return
        for i, rec in enumerate(records, 1):
            app = rec.app_name or ""
            if rec.app_version and rec.app_version != "unknown":
                app = f"{app} {rec.app_version}".strip()
            replay = f"  (replay of {rec.replayed_from})" if rec.replayed_from else ""
            started = rec.start_time[:16].replace("T", " ")
            if (duration := rec.duration) is not None:
                started += f" ({format_duration(duration)})"
            print(
                f"[{i}] {rec.run_id}  {rec.status:<11} {started}  "
                f"{rec.user or '?'}@{rec.host}  {app}{replay}"
            )
            if summary := bids_prov.source_summary(rec.source):
                print(f"    code: {summary}")
            print(f"    cd {shlex.quote(rec.cwd)}")
            print(f"    {rec.command}\n")

    @staticmethod
    def _exit(msg: str) -> NoReturn:
        print(msg, file=sys.stderr)
        sys.exit(1)
