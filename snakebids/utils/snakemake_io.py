"""File globbing functions based on snakemake.io library."""
from __future__ import annotations

import collections
import os
import re
from itertools import chain
from pathlib import Path
from typing import Sequence

from snakebids.types import ZipList
from snakebids.utils.containers import MultiSelectDict


def regex(filepattern: str) -> str:
    """Build Snakebids regex based on the given file pattern."""
    regex_list: list[str] = []
    last = 0
    wildcards: set[str] = set()
    for match in _wildcard_regex.finditer(filepattern):
        regex_list.append(re.escape(filepattern[last : match.start()]))
        wildcard = match.group("name")
        if wildcard in wildcards:
            if match.group("constraint"):
                msg = (
                    "Constraint regex must be defined only in the first "
                    "occurence of the wildcard in a string."
                )
                raise ValueError(msg)
            regex_list.append(f"(?P={wildcard})")
        else:
            wildcards.add(wildcard)
            constraint = (
                match.group("constraint") if match.group("constraint") else ".+"
            )
            regex_list.append(f"(?P<{wildcard}>{constraint})")

        last = match.end()
    regex_list.append(re.escape(filepattern[last:]))
    regex_list.append("$")  # ensure that the match spans the whole file
    return "".join(regex_list)


_wildcard_regex = re.compile(
    r"""
    \{
        (?=(   # This lookahead assertion emulates an 'atomic group'
               # which is required for performance
            \s*(?P<name>\w+)                    # wildcard name
            (\s*,\s*
                (?P<constraint>                 # an optional constraint
                    ([^{}]+ | \{\d+(,\d+)?\})*  # allow curly braces to nest
                )                               # ...  as in '{w,a{3,5}}'
            )?\s*
        ))\1
    \}
    """,
    re.VERBOSE,
)


def glob_wildcards(
    pattern: str | Path,
    files: Sequence[str | Path] | None = None,
    followlinks: bool = False,
) -> ZipList:
    """Glob the values of wildcards by matching a pattern to the filesystem.

    Returns a zip_list of field names with matched wildcard values.

    Parameters
    ----------
    pattern
        Path including wildcards to glob on the filesystem.
    files
        Files from which to glob wildcards. If None (default), the directory
        corresponding to the first wildcard in the pattern is walked, and
        wildcards are globbed from all files.
    followlinks
        Whether to follow links when globbing wildcards.
    """
    pattern = os.path.normpath(pattern)
    first_wildcard = re.search("{[^{]", pattern)
    dirname = (
        os.path.dirname(pattern[: first_wildcard.start()])
        if first_wildcard
        else os.path.dirname(pattern)
    )
    if not dirname:
        dirname = Path(".")

    names = [match.group("name") for match in _wildcard_regex.finditer(pattern)]

    # remove duplicates while preserving ordering
    names = list(dict.fromkeys(names))

    wildcards: dict[str, list[str]] = collections.defaultdict(list)

    re_pattern = re.compile(regex(pattern))

    file_iter = (
        (
            Path(dirpath, f)
            for dirpath, dirnames, filenames in os.walk(
                dirname, followlinks=followlinks
            )
            for f in chain(filenames, dirnames)
        )
        if files is None
        else iter(files)
    )

    for f in file_iter:
        if match := re.match(re_pattern, str(f)):
            for name, value in match.groupdict().items():
                wildcards[name].append(value)
    return MultiSelectDict(wildcards)


def update_wildcard_constraints(
    pattern: str,
    wildcard_constraints: dict[str, str],
    global_wildcard_constraints: dict[str, str],
) -> str:
    """Update wildcard constraints.

    Parameters
    ----------
    pattern : str
        Pattern on which to update constraints.
    wildcard_constraints : dict
        Dictionary of wildcard:constraint key-value pairs.
    global_wildcard_constraints : dict
        Dictionary of wildcard:constraint key-value pairs.
    """

    def replace_constraint(match: re.Match[str]):
        name = match.group("name")
        constraint = match.group("constraint")
        newconstraint = wildcard_constraints.get(
            name, global_wildcard_constraints.get(name)
        )
        if name in examined_names:
            return match.group(0)
        examined_names.add(name)
        # Don't override if constraint already set
        if constraint is not None:
            return match.group(0)
        # Only update if a new constraint has actually been set
        if newconstraint is not None:
            return f"{{{name},{newconstraint}}}"
        return match.group(0)

    examined_names: set[str] = set()
    return _wildcard_regex.sub(replace_constraint, pattern)
