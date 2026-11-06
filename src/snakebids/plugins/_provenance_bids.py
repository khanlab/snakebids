"""BIDS-facing views of snakebids run records.

The per-run records under ``code/snakebids/runs/`` are the single source of truth.
``dataset_description.json`` is *derived* from them and regenerated on every run.
It is kept as a compact summary: ``GeneratedBy`` gets one entry per distinct app
version (legacy "pipeline object" form, which current validators understand), not
one per run. Entries not written by snakebids are preserved, as are all other keys.

Also contains helpers to detect where the app's code came from (git checkout,
VCS/editable pip install) and whether it is running in a container.
"""

from __future__ import annotations

import importlib.metadata as impm
import json
import logging
import os
import re
import subprocess as sp
from collections.abc import Iterable, Sequence
from pathlib import Path
from typing import TYPE_CHECKING, Any, cast
from urllib.parse import unquote, urlparse

if TYPE_CHECKING:
    from snakebids.plugins.provenance import RunRecord

logger = logging.getLogger(__name__)

DEFAULT_BIDS_VERSION = "1.10.0"

#: Marker identifying GeneratedBy entries managed by snakebids
_MARKER = "Recorded by snakebids: "
_MARKER_RE = re.compile(r"^(Recorded by snakebids: |snakebids run )")


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
    return bool(_MARKER_RE.match(str(description)))


def _date(timestamp: str) -> str:
    return timestamp[:10]


def generated_by_summary(
    records: Sequence[RunRecord],
    name: str,
    provenance_dir: str,
    template: dict[str, Any] | None = None,
) -> list[dict[str, Any]]:
    """Summarize runs as ``GeneratedBy`` pipeline objects.

    One entry per distinct app version (not per run), so the list only grows when
    the software changes. Code URL and container are taken from the latest run of
    each version. The full details live in ``provenance_dir``.
    """
    groups: dict[str | None, list[RunRecord]] = {}
    for rec in records:
        version = rec.app_version if rec.app_version != "unknown" else None
        groups.setdefault(version, []).append(rec)

    entries: list[dict[str, Any]] = []
    for version, recs in groups.items():
        latest = recs[-1]
        first, last = _date(recs[0].start_time), _date(latest.start_time)
        n_failed = sum(r.status in {"failed", "interrupted"} for r in recs)
        description = (
            f"{_MARKER}{len(recs)} run{'s' if len(recs) > 1 else ''}"
            + (f" ({n_failed} failed)" if n_failed else "")
            + (f", {first}" if first == last else f", {first} to {last}")
        )
        if summary := source_summary(latest.source):
            description += f"; latest: {summary}"
        description += (
            f". Commands and configs are in {provenance_dir}/ "
            f"(list with: {name} --history <output_dir>)"
        )
        entry: dict[str, Any] = dict(template or {})
        entry.update({"Name": name, "Description": description})
        if version:
            entry["Version"] = version
        else:
            entry.pop("Version", None)
        if url := code_url(latest.source):
            entry["CodeURL"] = url
        if latest.container:
            entry["Container"] = latest.container
        entries.append(entry)
    return entries


def dataset_links(
    records: Iterable[RunRecord], existing: dict[str, str]
) -> dict[str, str]:
    """Add the input datasets of ``records`` to ``DatasetLinks``."""
    links = dict(existing)
    known = set(links.values())
    for record in records:
        for key, path in record.inputs.items():
            uri = Path(path).as_uri()
            if uri in known:
                continue
            alias, i = key, 1
            while alias in links:
                i += 1
                alias = f"{key}{i}"
            links[alias] = uri
            known.add(uri)
    return links


def update_dataset_description(
    output_dir: Path,
    records: Sequence[RunRecord],
    name: str,
    template: dict[str, Any] | None,
    provenance_dir: str = "code/snakebids",
) -> Path:
    """Create or update ``dataset_description.json`` with a summary of the runs.

    Parameters
    ----------
    template
        Contents of the app's ``pipeline_description.json``, used as the starting
        point when the dataset_description doesn't exist yet. A ``GeneratedBy`` entry
        in the template matching ``name`` supplies default keys (e.g. ``Author``)
        for the summary entries.
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

    if existing and all(isinstance(e, str) for e in existing):
        # GeneratedBy uses (draft BEP028) activity identifiers; mixing in pipeline
        # objects would make it invalid, so leave it alone.
        logger.warning(
            "%s uses identifiers in GeneratedBy; not adding %s runs to it", path, name
        )
    else:
        kept = [
            e
            for e in existing
            if not _is_snakebids_entry(e) and not (app_tmpl and e == app_tmpl)
        ]
        desc["GeneratedBy"] = kept + generated_by_summary(
            records, name, provenance_dir, app_tmpl
        )

    links = desc.get("DatasetLinks", {})
    links = dataset_links(
        records,
        {str(k): str(v) for k, v in cast("dict[Any, Any]", links).items()}
        if isinstance(links, dict)
        else {},
    )
    if links:
        desc["DatasetLinks"] = links
    write_json(path, desc)
    return path
