"""Tests for snakebids.utils.snakemake_templates module."""

from __future__ import annotations

import re

import pytest

from snakebids.utils.snakemake_templates import SnakemakeWildcards


def find_non_empty_match(pattern: re.Pattern[str], text: str) -> str:
    """Find the first non-empty match in text.

    Since wildcard patterns are optional (ending with ?), they always match,
    but may match empty string. This helper finds the first non-empty match.

    Returns
    -------
    str
        The matched string, or empty string if no non-empty match found.
    """
    for match in pattern.finditer(text):
        if match.group():
            return match.group()
    return ""


def do_match(pattern: str, text: str) -> str | None:
    match = re.match(f"(?:{pattern})$", text)
    if match is not None:
        return match.group(1)
    return None


def test_directory_formats_label_correctly():
    """Test that directory wildcard label is formatted correctly."""
    wc = SnakemakeWildcards("run")
    name, _ = wc.directory.split(",", 1)
    assert name == "_run_d_"


@pytest.mark.parametrize(
    ("before", "text", "after", "expected"),
    [
        ("", "", "", ""),
        ("", "more", "more", ""),
        ("prefix", "prefix", "", ""),
        ("prefix", "prefixmore", "more", ""),
        ("", "run-01/", "", "run-01/"),
        ("", "nur-01/", "", None),
        ("", "run-/", "", None),
        ("", "run01/", "", None),
        ("", "run-_/", "", None),
        ("", "run--/", "", None),
        ("", "run-\n/", "", None),
        ("", "run-01//", "", None),
        ("", "run-01", "", None),
        ("", "/run-01/", "", None),
        ("/", "/run-01/", "", "run-01/"),
        ("", "prefixrun-01/", "", None),
        ("prefix", "prefixrun-01/", "", None),
        ("prefix", "prefix/run-01/", "", None),
        ("prefix/", "prefix/run-01/", "", "run-01/"),
        ("prefix/", "prefix/run-01/more", "more", "run-01/"),
        ("", "run-01_test/", "", None),
        ("", "run-01-02/", "", None),
    ],
)
def test_directory_wildcard_matching(before: str, after: str, text: str, expected: str):
    """Test directory wildcard matching behavior."""
    wc = SnakemakeWildcards("run")
    _, constraint = wc.directory.split(",", 1)
    assert do_match(f"{before}({constraint}){after}", text) == expected


@pytest.mark.parametrize(
    ("before", "text", "after", "expected"),
    [
        ("prefix", "prefix/", "", "/"),
        ("", "prefix/", "", None),
        ("", "", "", ""),
        ("", "more", "more", ""),
        ("prefix", "prefixmore", "more", None),
        ("prefix", "prefix/more", "more", "/"),
        ("prefix/", "prefix/more", "more", ""),
        ("prefix", "prefix//more", "more", None),
        ("", "/more", "more", None),
        ("/", "/more", "more", ""),
        ("", "/", "", None),
        ("/", "/", "", ""),
        ("/", "//", "", "/"),
    ],
)
def test_d_wildcard_matching(before: str, after: str, text: str, expected: bool):
    """Test __d__ wildcard matching behavior."""
    _, constraint = SnakemakeWildcards.d.split(",", 1)
    assert do_match(f"{before}({constraint}){after}", text) == expected


@pytest.mark.parametrize(
    ("before", "text", "after", "expected"),
    [
        ("", "", "", ""),
        ("", "more", "more", ""),
        ("prefix", "prefix", "", None),
        ("prefix", "prefixmore", "more", None),
        ("prefix", "prefix_more", "more", "_"),
        ("prefix", "prefix_.more", "more", None),
        ("prefix", "prefix/more", "more", None),
        ("prefix/", "prefix/more", "more", ""),
        ("prefix/", "prefix/_more", "more", None),
        ("prefix", "prefix/_more", "more", None),
        ("", "_more", "more", None),
        ("", "_m", "", None),
        ("", "_.", "\\.", None),
        ("", "_", "", None),
    ],
)
def test_underscore_wildcard_matching(
    before: str, after: str, text: str, expected: bool
):
    _, constraint = SnakemakeWildcards.underscore.split(",", 1)
    assert do_match(f"{before}({constraint}){after}", text) == expected


@pytest.mark.parametrize(
    ("before", "text", "after", "expected"),
    [
        ("", "", "", ""),
        ("", "more", "more", ""),
        ("prefix", "prefix", "", ""),
        ("prefix", "prefixmore", "more", ""),
        ("", "anat01/", "/", "anat01"),
        ("", "anat/", "", None),
        ("", "anat", "", None),
        ("", "anat_/", "/", None),
        ("", "anat-01/", "/", None),
        ("", "an/a/", "/", None),
        ("", "ana\n/", "/", None),
        ("", "/ana-01/", "/", None),
        ("/", "/anat01/", "", None),
        ("/", "/anat01/", "/", "anat01"),
        ("prefix", "prefixanat01/", "/", None),
        ("prefix", "prefix/anat/", "/", None),
        ("prefix/", "prefix/anat/", "", None),
        ("prefix/", "prefix/anat/", "/", "anat"),
        ("prefix/", "prefix/anat01/more", "/more", "anat01"),
    ],
)
def test_datatype_wildcard_matching(before: str, after: str, text: str, expected: bool):
    _, constraint = SnakemakeWildcards.datatype.split(",", 1)
    assert do_match(f"{before}({constraint}){after}", text) == expected


@pytest.mark.parametrize(
    ("before", "text", "after", "expected"),
    [
        ("", "", "", ""),
        ("", "more", "more", ""),
        ("prefix", "prefix", "", ""),
        ("prefix", "prefixmore", "more", ""),
        ("", "run-0", "0", "run-"),
        ("", "nur-0", "0", None),
        ("", "run-", "", None),
        ("", "run-0", "", None),
        ("", "run0", "0", None),
        ("", "run", "", None),
        ("", "run-_", "_", None),
        ("", "run--", "-", None),
        ("", "run-\n", "\n", None),
        ("", "run-/", "/", None),
        ("", "/run-0", "0", None),
        ("/", "/run-0", "0", "run-"),
        ("_", "_run-0", "0", None),
        ("", "_run-0", "0", None),
        ("prefix_", "prefix_run-0", "0", None),
        ("prefix", "prefix_run-0", "0", "_run-"),
        ("prefix", "prefix/run-0", "0", None),
        ("prefix/", "prefix/run-0", "0", "run-"),
        ("prefix/", "prefix/_run-0", "0", None),
        ("prefix/_", "prefix/_run-0", "0", None),
        ("prefix/", "prefix/run-0/more", "0/more", "run-"),
    ],
)
def test_dummy_wildcard_matching(before: str, after: str, text: str, expected: bool):
    _, constraint = SnakemakeWildcards("run").dummy.split(",", 1)
    assert do_match(f"{before}({constraint}){after}", text) == expected


@pytest.mark.parametrize(
    ("before", "text", "after", "expected"),
    [
        ("", "", "", ""),
        ("", "more", "more", ""),
        ("prefix", "prefix", "", ""),
        ("prefix", "prefixmore", "more", ""),
        ("run-", "run-value0", "", "value0"),
        ("", "run-value0", "", None),
        ("run-", "nur-value0", "", None),
        ("run-", "run-", "", ""),
        ("run", "runvalue0", "", None),
        ("run-", "run-_", "", None),
        ("run-", "run--", "", None),
        ("run-", "run-\nfwaef", "", None),
        ("run-", "run-/", "", None),
        ("/run-", "/run-value0", "", "value0"),
        ("prefix/run-", "prefix/run-0", "", "0"),
        ("prefixrun-", "prefixrun-0", "", "0"),
        ("run-", "run-00", "0", None),
        ("run-", "run-0more", "more", None),
        ("run-", "run-0_more", "_more", "0"),
        ("run-", "run-0/more", "/more", "0"),
        ("run-", "run-0_more", "more", None),
        ("run-", "run-0/more", "more", None),
    ],
)
def test_ordinary_wildcard_matching(before: str, after: str, text: str, expected: bool):
    _, constraint = SnakemakeWildcards("run").variable.split(",", 1)
    assert do_match(f"{before}({constraint}){after}", text) == expected


@pytest.mark.parametrize(
    ("before", "text", "after", "expected"),
    [
        ("", "", "", ""),
        ("", "more", "more", ""),
        ("prefix", "prefix", "", ""),
        ("prefix", "prefixmore", "more", ""),
        ("", "suffix.ext", "", "suffix.ext"),
        ("", "suffix-ext", "", None),
        ("", "suffix_ext", "", None),
        ("", "suffix/ext", "", None),
        ("", "_suffix.ext", "", None),
        ("_", "_suffix.ext", "", "suffix.ext"),
        ("/_", "/_suffix.ext", "", "suffix.ext"),
        ("", "/suffix.ext", "", None),
        ("/", "/suffix.ext", "", "suffix.ext"),
        ("/", "/suffix.extmore", "more", "suffix.ext"),
        ("/", "/suffix.ext_more", "_more", "suffix.ext"),
        ("/", "/suffix.ext/more", "/more", "suffix.ext"),
        ("/", "/suffix.ext-more", "-more", "suffix.ext"),
    ],
)
def test_suffix_wildcard_matching(before: str, after: str, text: str, expected: bool):
    _, constraint = SnakemakeWildcards.suffix.split(",", 1)
    assert do_match(f"{before}({constraint}){after}", text) == expected


@pytest.mark.parametrize(
    ("before", "text", "after", "expected"),
    [
        ("", "", "", ""),
        ("", "more", "more", ""),
        ("prefix", "prefix", "", ""),
        ("prefix", "prefixmore", "more", ""),
        ("", ".ext", "", ".ext"),
        ("", ".ext-more", "", None),
        ("", ".ext_more", "", None),
        ("", ".ext/more", "", None),
        ("", "before.ext", "", None),
        ("before", "before.ext", "", ".ext"),
        ("/_", "/_.ext", "", ".ext"),
        ("/_", "/_.extmore", "more", None),
        ("/_", "/_.ext/more", "/more", None),
        ("/_", "/_.ext_more", "_more", None),
    ],
)
def test_extension_wildcard_matching(
    before: str, after: str, text: str, expected: bool
):
    _, constraint = SnakemakeWildcards.extension.split(",", 1)
    assert do_match(f"{before}({constraint}){after}", text) == expected


@pytest.mark.parametrize(
    ("tag", "expected_name"),
    [
        ("sub", "subject"),
        ("ses", "session"),
        ("run", "run"),
    ],
)
def test_tag_replacement(tag: str, expected_name: str):
    """Test tag to wildcard name conversion."""
    wc = SnakemakeWildcards(tag)
    assert expected_name in wc.variable
    assert expected_name in wc.dummy
    assert expected_name in wc.directory
