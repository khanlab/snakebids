"""BIDS-facing views of snakebids run records.

The per-run records under ``code/snakebids/runs/`` are the single source of truth.
Everything here is *derived* from them and regenerated on every run:

- ``dataset_description.json``: ``GeneratedBy`` gets one entry per recorded run
  (legacy "pipeline object" form, which current validators understand). Entries not
  written by snakebids are preserved, as are all other keys.
- ``prov/prov-<app>_{act,soft,env,ent}.json`` and ``prov/provenance.tsv``: loosely
  following the draft BEP028 provenance spec (W3C PROV-inspired). This is draft and
  likely to change; it is written for inspiration/forward-compatibility only.

Also contains helpers to detect where the app's code came from (git checkout,
VCS/editable pip install) and whether it is running in a container.
"""

from __future__ import annotations

import csv
import hashlib
import importlib.metadata as impm
import json
import logging
import os
import re
import subprocess as sp
from collections.abc import Iterable, Sequence
from datetime import datetime
from pathlib import Path
from typing import TYPE_CHECKING, Any, cast
from urllib.parse import unquote, urlparse

if TYPE_CHECKING:
    from snakebids.plugins.provenance import RunRecord

logger = logging.getLogger(__name__)

DEFAULT_BIDS_VERSION = "1.10.0"

#: Marker identifying GeneratedBy entries managed by snakebids
_RUN_MARKER = "snakebids run "
_RUN_MARKER_RE = re.compile(rf"^{_RUN_MARKER}(\S+)")


# --------------------------------------------------------------------------------------
# Code source and container detection
# --------------------------------------------------------------------------------------


def _git(args: Sequence[str], cwd: Path) -> str | None:
    try:
        proc = sp.run(
            ["git", *args],
            cwd=cwd,
            capture_output=True,
            text=True,
            timeout=10,
            check=False,
        )
    except (OSError, sp.SubprocessError):
        return None
    if proc.returncode != 0:
        return None
    return proc.stdout.strip()


def sanitize_url(url: str) -> str:
    """Strip credentials and normalize ssh remotes to https where unambiguous."""
    # scp-like ssh syntax: git@github.com:org/repo.git
    if m := re.fullmatch(r"[\w.-]+@([\w.-]+):(?!//)(.+?)(?:\.git)?/?", url):
        return f"https://{m[1]}/{m[2]}"
    url = re.sub(r"^(\w[\w+.-]*://)[^/@]+@", r"\1", url)  # drop user:token@
    if url.startswith("ssh://"):
        url = "https://" + url.removeprefix("ssh://")
    return re.sub(r"\.git/?$", "", url)


def git_info(path: Path) -> dict[str, Any] | None:
    """Describe the git checkout containing ``path``, if its files are tracked.

    Returns None if ``path`` is not inside a git work tree, or if it is not tracked
    by that repository (e.g. the app is pip-installed into a venv that happens to
    live inside some unrelated clone).
    """
    path = Path(path)
    if not path.is_dir():
        return None
    if _git(["rev-parse", "--is-inside-work-tree"], path) != "true":
        return None
    if not _git(["ls-files", "--", "."], path):
        return None
    commit = _git(["rev-parse", "HEAD"], path)
    if commit is None:
        return None
    info: dict[str, Any] = {
        "commit": commit,
        "toplevel": _git(["rev-parse", "--show-toplevel"], path),
        "describe": _git(["describe", "--tags", "--always", "--dirty"], path),
        "branch": _git(["rev-parse", "--abbrev-ref", "HEAD"], path),
        "dirty": bool(_git(["status", "--porcelain", "--untracked-files=no"], path)),
    }
    if remote := _git(["config", "--get", "remote.origin.url"], path):
        info["remote"] = sanitize_url(remote)
    return info


_MAX_DIFF_BYTES = 10 * 1024 * 1024


def write_git_diff(repo: Path, dest: Path) -> bool:
    """Save uncommitted changes to tracked files, so a dirty run can be reproduced.

    Apply with ``git apply <dest>`` on top of the recorded commit.
    """
    diff = _git(["diff", "HEAD", "--binary"], repo)
    if not diff:
        return False
    if len(diff) > _MAX_DIFF_BYTES:
        logger.warning("Uncommitted changes too large to record (%d bytes)", len(diff))
        return False
    dest.write_text(diff + "\n")
    return True


def _direct_url(distribution: str) -> dict[str, Any] | None:
    """Read the PEP 610 ``direct_url.json`` of an installed distribution."""
    try:
        text = impm.distribution(distribution).read_text("direct_url.json")
    except (impm.PackageNotFoundError, OSError):
        return None
    if not text:
        return None
    try:
        return json.loads(text)
    except json.JSONDecodeError:
        return None


def guess_distribution(app_dir: Path | None) -> str | None:
    """Find the installed distribution providing the package at ``app_dir``."""
    if app_dir is None:
        return None
    dists = impm.packages_distributions().get(app_dir.name, [])
    return dists[0] if len(dists) == 1 else None


def code_source(
    app_dir: Path | None, distribution: str | None
) -> dict[str, Any] | None:
    """Work out where the running app code came from.

    In order of preference:

    1. ``app_dir`` is a tracked part of a git checkout (``git clone`` + run)
    2. the distribution is an editable install of a git checkout
    3. the distribution was pip-installed from a VCS URL (commit recorded by pip)
    4. the distribution was installed from an index/wheel (nothing extra to record)
    """
    if app_dir is not None and (info := git_info(app_dir)):
        return {"type": "git", **info}
    if distribution is None or (direct := _direct_url(distribution)) is None:
        return None
    url = str(direct.get("url", ""))
    if direct.get("dir_info", {}).get("editable") and url.startswith("file://"):
        src_path = Path(unquote(urlparse(url).path))
        if info := git_info(src_path):
            return {"type": "git", "editable": True, **info}
        return {"type": "editable", "path": str(src_path)}
    if vcs := direct.get("vcs_info"):
        return {
            "type": "vcs",
            "remote": sanitize_url(url),
            "commit": vcs.get("commit_id"),
            "requested_revision": vcs.get("requested_revision"),
        }
    return {"type": "url", "url": sanitize_url(url)}


def code_url(source: dict[str, Any] | None) -> str | None:
    """Best-effort browsable URL for the exact code that was run."""
    if not source or not (remote := source.get("remote")):
        return None
    commit = source.get("commit")
    if (
        commit
        and remote.startswith("https://")
        and re.match(r"https://(github\.com|gitlab\.[\w.]+|bitbucket\.org)/", remote)
    ):
        return f"{remote}/tree/{commit}"
    return remote


def source_summary(source: dict[str, Any] | None) -> str | None:
    """Short human-readable description of the code source."""
    if not source:
        return None
    if source.get("type") in {"git", "vcs"} and source.get("commit"):
        ref = source.get("describe") or source["commit"][:12]
        dirty = " (uncommitted changes, see app.diff)" if source.get("dirty") else ""
        where = source.get("remote") or source.get("path") or ""
        return f"git {ref}{dirty} {where}".strip()
    if source.get("type") == "editable":
        return f"editable install at {source.get('path')}"
    return source.get("url")


def detect_container() -> dict[str, Any] | None:
    """Detect apptainer/singularity/docker, returning a GeneratedBy Container object."""
    for var, kind in (
        ("APPTAINER_CONTAINER", "apptainer"),
        ("SINGULARITY_CONTAINER", "singularity"),
    ):
        if image := os.environ.get(var):
            name = os.environ.get(var.replace("_CONTAINER", "_NAME"))
            container = {"Type": kind, "Tag": name or Path(image).name}
            if "://" in image:
                container["URI"] = image
            return container
    if Path("/.dockerenv").exists():
        return {"Type": "docker"}
    return None


# --------------------------------------------------------------------------------------
# dataset_description.json
# --------------------------------------------------------------------------------------


def json_objects(value: Any) -> list[dict[str, Any]]:
    """Filter a json value down to its list of objects."""
    if not isinstance(value, list):
        return []
    return [
        cast("dict[str, Any]", e)
        for e in cast("list[Any]", value)
        if isinstance(e, dict)
    ]


def read_json(path: Path) -> dict[str, Any] | None:
    try:
        data = json.loads(path.read_text())
    except FileNotFoundError:
        return None
    except (OSError, json.JSONDecodeError) as err:
        logger.warning("Could not read %s: %s", path, err)
        return None
    return cast("dict[str, Any]", data) if isinstance(data, dict) else None


def write_json(path: Path, data: Any) -> None:
    """Atomically write json to ``path``."""
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(f".{path.name}.tmp")
    tmp.write_text(json.dumps(data, indent=2) + "\n")
    tmp.replace(path)


def _is_snakebids_entry(entry: Any) -> bool:
    if not isinstance(entry, dict):
        return False
    description = cast("dict[str, Any]", entry).get("Description", "")
    return bool(_RUN_MARKER_RE.match(str(description)))


def generated_by_entry(
    record: RunRecord, name: str, template: dict[str, Any] | None = None
) -> dict[str, Any]:
    """Legacy-style ``GeneratedBy`` pipeline object for one run."""
    description = f"{_RUN_MARKER}{record.run_id} [{record.status}]"
    if record.replayed_from:
        description += f" (replay of {record.replayed_from})"
    if summary := source_summary(record.source):
        description += f"; {summary}"
    description += f": {record.command}"
    entry: dict[str, Any] = dict(template or {})
    entry.update({"Name": name, "Description": description})
    if record.app_version and record.app_version != "unknown":
        entry["Version"] = record.app_version
    if url := code_url(record.source):
        entry["CodeURL"] = url
    if record.container:
        entry["Container"] = record.container
    return entry


def existing_dataset_links(
    output_dir: Path, template: dict[str, Any] | None
) -> dict[str, str]:
    """DatasetLinks of the current dataset_description, or else of the template."""
    desc = read_json(output_dir / "dataset_description.json")
    if desc is None:
        desc = template or {}
    links = desc.get("DatasetLinks", {})
    if not isinstance(links, dict):
        return {}
    return {str(k): str(v) for k, v in cast("dict[Any, Any]", links).items()}


def dataset_links(
    records: Iterable[RunRecord], existing: dict[str, str]
) -> tuple[dict[str, str], dict[str, str]]:
    """Assign DatasetLinks aliases to input datasets.

    Returns updated links and a mapping of input path -> alias.
    """
    links = dict(existing)
    by_uri = {uri: alias for alias, uri in links.items()}
    aliases: dict[str, str] = {}
    for record in records:
        for key, path in record.inputs.items():
            uri = Path(path).as_uri()
            if uri in by_uri:
                aliases[path] = by_uri[uri]
                continue
            alias, i = key, 1
            while alias in links:
                i += 1
                alias = f"{key}{i}"
            links[alias] = uri
            by_uri[uri] = alias
            aliases[path] = alias
    return links, aliases


def update_dataset_description(
    output_dir: Path,
    records: Sequence[RunRecord],
    name: str,
    template: dict[str, Any] | None,
    activity_ids: dict[str, str] | None = None,
) -> Path:
    """Create or update ``dataset_description.json`` with the recorded runs.

    Parameters
    ----------
    template
        Contents of the app's ``pipeline_description.json``, used as the starting
        point when the dataset_description doesn't exist yet. A ``GeneratedBy`` entry
        in the template matching ``name`` supplies default keys (e.g. ``Author``)
        for each run entry.
    activity_ids
        run_id -> BEP028 activity Id. Used if the existing file already uses the
        BEP028 identifier form of ``GeneratedBy`` (a list of strings).
    """
    path = output_dir / "dataset_description.json"
    template = template or {}
    tmpl_entries = json_objects(template.get("GeneratedBy"))
    app_tmpl = next((e for e in tmpl_entries if e.get("Name") == name), None)

    desc = read_json(path)
    if desc is None:
        desc = {k: v for k, v in template.items() if k != "GeneratedBy"}
        desc.setdefault("Name", f"Dataset generated by {name}")
        desc.setdefault("BIDSVersion", DEFAULT_BIDS_VERSION)
        desc.setdefault("DatasetType", "derivative")
        existing: list[Any] = [e for e in tmpl_entries if e is not app_tmpl]
    else:
        raw = desc.get("GeneratedBy", [])
        existing = [raw] if isinstance(raw, str | dict) else list(raw)

    if existing and all(isinstance(e, str) for e in existing) and activity_ids:
        # BEP028 identifier form: keep foreign ids, (re)add ours in order
        ours = set(activity_ids.values())
        desc["GeneratedBy"] = [e for e in existing if e not in ours] + [
            activity_ids[r.run_id] for r in records if r.run_id in activity_ids
        ]
    else:
        kept = [
            e
            for e in existing
            if not _is_snakebids_entry(e) and not (app_tmpl and e == app_tmpl)
        ]
        desc["GeneratedBy"] = kept + [
            generated_by_entry(r, name, app_tmpl) for r in records
        ]

    links, _ = dataset_links(records, desc.get("DatasetLinks", {}))
    if links:
        desc["DatasetLinks"] = links
    write_json(path, desc)
    return path


# --------------------------------------------------------------------------------------
# BEP028 (draft) provenance files
# --------------------------------------------------------------------------------------


def _uid(obj: Any) -> str:
    """Content hash, so ids differ exactly when attributes differ (as BEP028 asks)."""
    blob = json.dumps(obj, sort_keys=True, default=str).encode()
    return hashlib.sha256(blob).hexdigest()[:8]


def _with_id(prefix: str, obj: dict[str, Any], extra: Any = None) -> dict[str, Any]:
    """Prepend an Id derived from the object's content (plus ``extra`` hidden data)."""
    return {"Id": f"bids::prov#{prefix}-{_uid([obj, extra])}", **obj}


def _seconds(timestamp: str) -> str:
    try:
        return datetime.fromisoformat(timestamp).isoformat(timespec="seconds")
    except ValueError:
        return timestamp


def prov_label(name: str) -> str:
    """Sanitize an app name for use as a BIDS ``prov-<label>``."""
    return re.sub(r"[^a-zA-Z0-9]", "", name) or "snakebids"


def build_prov(
    records: Sequence[RunRecord],
    name: str,
    dataset_links: dict[str, str],
    provenance_dir: str,
) -> tuple[dict[str, dict[str, Any]], dict[str, str]]:
    """Build BEP028-style provenance documents from run records.

    Returns a mapping of file suffix (``act``, ``soft``, ``env``, ``ent``) to json
    document, and a mapping of run_id -> activity Id.
    """
    label = prov_label(name)
    software: dict[str, dict[str, Any]] = {}
    environments: dict[str, dict[str, Any]] = {}
    datasets: dict[str, dict[str, Any]] = {}
    activities: list[dict[str, Any]] = []
    activity_ids: dict[str, str] = {}
    by_uri = {uri: alias for alias, uri in dataset_links.items()}

    def add(store: dict[str, dict[str, Any]], obj: dict[str, Any]) -> str:
        store.setdefault(obj["Id"], obj)
        return obj["Id"]

    for rec in records:
        app_sw: dict[str, Any] = {"Label": name, "Version": rec.app_version}
        if url := code_url(rec.source):
            app_sw["AlternativeIdentifier"] = [url]
        # the exact code (commit + dirty state) distinguishes software objects
        src = rec.source or {}
        app_id = add(
            software,
            _with_id(label, app_sw, [src.get("commit"), src.get("dirty")]),
        )
        sw_ids = [app_id]
        for dep, version in (
            ("snakebids", rec.snakebids_version),
            ("snakemake", rec.snakemake_version),
        ):
            if version:
                sw_ids.append(
                    add(
                        software,
                        _with_id(
                            dep,
                            {
                                "Label": dep,
                                "Version": version,
                                "ActedOnBehalfOf": [app_id],
                            },
                        ),
                    )
                )

        env: dict[str, Any] = {
            "Label": (rec.container or {}).get("Tag") or rec.host,
            "OperatingSystem": rec.operating_system,
            "Dependencies": {"python": rec.python_version},
        }
        if uri := (rec.container or {}).get("URI"):
            env["AlternativeIdentifier"] = [uri]
        env = {k: v for k, v in env.items() if v}
        used = [add(environments, _with_id("env", env))]

        for path in rec.inputs.values():
            alias = by_uri.get(Path(path).as_uri())
            if alias is None:
                continue
            ds_id = f"bids:{alias}:."
            datasets.setdefault(ds_id, {"Id": ds_id, "Label": _dataset_name(path)})
            used.append(ds_id)
        used.append(f"bids::{provenance_dir}/runs/{rec.run_id}/config.yaml")

        act_id = f"bids::prov#{label}-{rec.run_id}"
        activity_ids[rec.run_id] = act_id
        act: dict[str, Any] = {
            "Id": act_id,
            "Label": f"{name} run {rec.run_id}",
            "Command": rec.command,
            "Description": (
                f"status: {rec.status}; working directory: {rec.cwd}; "
                f"host: {rec.host}" + (f"; user: {rec.user}" if rec.user else "")
            ),
            "AssociatedWith": sw_ids,
            "Used": used,
            "StartedAtTime": _seconds(rec.start_time),
        }
        if rec.end_time:
            act["EndedAtTime"] = _seconds(rec.end_time)
        activities.append(act)

    docs: dict[str, dict[str, Any]] = {
        "act": {"Activities": activities},
        "soft": {"Software": list(software.values())},
        "env": {"Environments": list(environments.values())},
    }
    if datasets:
        docs["ent"] = {"Datasets": list(datasets.values())}
    return docs, activity_ids


def _dataset_name(path: str) -> str:
    desc = read_json(Path(path) / "dataset_description.json")
    if desc and isinstance(desc.get("Name"), str) and desc["Name"]:
        return desc["Name"]
    return Path(path).name


def write_prov(
    output_dir: Path,
    docs: dict[str, dict[str, Any]],
    name: str,
    prov_dir: str = "prov",
) -> None:
    """Write BEP028 files and register the label in ``prov/provenance.tsv``."""
    label = prov_label(name)
    root = output_dir / prov_dir
    for suffix in ("act", "soft", "env", "ent"):
        path = root / f"prov-{label}_{suffix}.json"
        if suffix in docs:
            write_json(path, docs[suffix])
        else:
            path.unlink(missing_ok=True)

    tsv = root / "provenance.tsv"
    prov_id = f"prov-{label}"
    rows: list[dict[str, str]] = []
    fields = ["provenance_id", "description"]
    if tsv.exists():
        with tsv.open(newline="") as f:
            reader = csv.DictReader(f, delimiter="\t")
            fields = list(reader.fieldnames or fields)
            rows = [r for r in reader if r.get("provenance_id") != prov_id]
    if "description" not in fields:
        fields.append("description")
    rows.append(
        {
            "provenance_id": prov_id,
            "description": f"Runs of {name}, recorded by snakebids",
        }
    )
    tmp = tsv.with_name(f".{tsv.name}.tmp")
    with tmp.open("w", newline="") as f:
        writer = csv.DictWriter(
            f, fieldnames=fields, delimiter="\t", restval="n/a", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)
    tmp.replace(tsv)
