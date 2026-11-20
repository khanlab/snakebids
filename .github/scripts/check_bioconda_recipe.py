# ruff: noqa: INP001
"""Check that the bioconda recipe's run requirements match pyproject.toml.

The bioconda autobump bot updates the recipe's version and checksum on each
release, but never its requirements, so they drift whenever pyproject.toml's
dependencies change. This compares the two and prints the `run:` block the
recipe should have.

Usage: check_bioconda_recipe.py PYPROJECT META_YAML [--summary FILE]

Exits 1 if the requirements differ.
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import tomllib
from packaging.requirements import Requirement
from packaging.specifiers import SpecifierSet
from packaging.utils import canonicalize_name

# pypi name -> conda-forge name, where they differ
CONDA_NAMES = {
    "typing-extensions": "typing_extensions",
    "ruamel-yaml": "ruamel.yaml",
    "universal-pathlib": "universal_pathlib",
}

# What to do with dependencies that have python_version markers. Noarch conda
# packages can't express these, so each one needs a decision:
#   "keep": use the spec that applies on the newest python
#   "drop": leave out of the recipe (e.g. wheel-availability hints for
#           transitive dependencies)
MARKER_POLICY = {
    "snakemake": "keep",
    "numpy": "drop",
    "pandas": "drop",
}
NEWEST_PYTHON = {"python_version": "3.99", "python_full_version": "3.99.0"}


def conda_name(name: str) -> str:
    """Map a pypi package name to its conda-forge name."""
    name = canonicalize_name(name)
    return CONDA_NAMES.get(name, name)


def expand(spec: SpecifierSet) -> frozenset[str]:
    """Normalize a specifier set, expanding ~= into its >=,< equivalent."""
    out = set()
    for s in spec:
        if s.operator == "~=":
            parts = s.version.split(".")
            upper = parts[:-1]
            upper[-1] = str(int(upper[-1]) + 1)
            out |= {f">={s.version}", f"<{'.'.join(upper)}"}
        else:
            out.add(f"{s.operator}{s.version}")
    return frozenset(out)


def from_pyproject(path: Path) -> tuple[str, dict[str, SpecifierSet]]:
    """Read python_min and the conda run requirements from pyproject.toml."""
    project = tomllib.loads(path.read_text())["project"]
    python_min = SpecifierSet(project["requires-python"])
    python_min = next(s.version for s in python_min if s.operator == ">=")

    deps: dict[str, SpecifierSet] = {}
    for raw in project["dependencies"]:
        req = Requirement(raw)
        name = conda_name(req.name)
        if req.marker is None:
            deps[name] = req.specifier
            continue
        policy = MARKER_POLICY.get(name)
        if policy is None:
            sys.exit(
                f"{raw!r} has an environment marker; add {name!r} to "
                "MARKER_POLICY in this script"
            )
        if policy == "keep" and req.marker.evaluate(NEWEST_PYTHON):
            deps[name] = req.specifier
    return python_min, deps


def from_recipe(path: Path) -> tuple[str | None, dict[str, SpecifierSet]]:
    """Read python_min and the run requirements from a recipe meta.yaml."""
    text = path.read_text()
    m = re.search(r'set python_min\s*=\s*"([^"]+)"', text)
    python_min = m.group(1) if m else None

    m = re.search(r"^  run:\n((?:    .*\n)+)", text, re.MULTILINE)
    if not m:
        sys.exit(f"no requirements/run block found in {path}")
    deps: dict[str, SpecifierSet] = {}
    for raw in m.group(1).splitlines():
        entry = raw.split("#", 1)[0].strip()
        if not entry.startswith("- "):
            continue
        name, _, spec = entry[2:].strip().partition(" ")
        if name == "python":
            continue
        deps[name] = SpecifierSet(spec.replace(" ", ""))
    return python_min, deps


def run_block(deps: dict[str, SpecifierSet]) -> str:
    """Format the requirements as a recipe `run:` block."""
    lines = ["  run:", "    - python >={{ python_min }}"]
    for name, spec in deps.items():
        spec_str = ",".join(sorted(expand(spec), key=lambda s: s.lstrip("<>=!~")))
        lines.append(f"    - {name} {spec_str}".rstrip())
    return "\n".join(lines)


def main() -> int:
    """Compare the requirements and report any differences."""
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("pyproject", type=Path)
    parser.add_argument("recipe", type=Path)
    parser.add_argument("--summary", type=Path, help="write a markdown report here")
    args = parser.parse_args()

    want_py, want = from_pyproject(args.pyproject)
    have_py, have = from_recipe(args.recipe)

    problems = []
    if have_py != want_py:
        problems.append(f"python_min: recipe has `{have_py}`, expected `{want_py}`")
    for name in sorted(want.keys() | have.keys()):
        if name not in have:
            problems.append(f"`{name}`: missing from recipe (expected `{want[name]}`)")
        elif name not in want:
            problems.append(f"`{name}`: in recipe but not in pyproject.toml")
        elif expand(have[name]) != expand(want[name]):
            problems.append(
                f"`{name}`: recipe has `{have[name]}`, expected `{want[name]}`"
            )

    if not problems:
        print("bioconda recipe requirements match pyproject.toml")
        return 0

    report = "\n".join(
        [
            "The bioconda recipe's requirements differ from `pyproject.toml`:",
            "",
            *(f"- {p}" for p in problems),
            "",
            (
                "Expected `requirements: run:` block "
                f'(with `{{% set python_min = "{want_py}" %}}`):'
            ),
            "",
            "```yaml",
            run_block(want),
            "```",
        ]
    )
    print(report)
    if args.summary:
        args.summary.write_text(report + "\n")
    return 1


if __name__ == "__main__":
    sys.exit(main())
