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


def do_match(pattern: str, text: str) -> bool:
    return re.match(f"({pattern})$", text) is not None


def test_directory_formats_label_correctly():
    """Test that directory wildcard label is formatted correctly."""
    wc = SnakemakeWildcards("run")
    name, _ = wc.directory.split(",", 1)
    assert name == "_run_d_"


@pytest.mark.parametrize(
    ("before", "text", "after", "expected"),
    [
        ("", "", "", True),
        ("", "more", "more", True),
        ("prefix", "prefix", "", True),
        ("prefix", "prefixmore", "more", True),
        ("", "run-01/", "", True),
        ("", "nur-01/", "", False),
        ("", "run-/", "", False),
        ("", "run01/", "", False),
        ("", "run-_/", "", False),
        ("", "run--/", "", False),
        ("", "run-\n/", "", False),
        ("", "run-01//", "", False),
        ("", "run-01", "", False),
        ("", "/run-01/", "", False),
        ("/", "/run-01/", "", True),
        ("", "prefixrun-01/", "", False),
        ("prefix", "prefixrun-01/", "", False),
        ("prefix", "prefix/run-01/", "", False),
        ("prefix/", "prefix/run-01/", "", True),
        ("prefix/", "prefix/run-01/more", "more", True),
        ("", "run-01_test/", "", False),
        ("", "run-01-02/", "", False),
    ],
)
def test_directory_wildcard_matching(
    before: str, after: str, text: str, expected: bool
):
    """Test directory wildcard matching behavior."""
    wc = SnakemakeWildcards("run")
    _, constraint = wc.directory.split(",", 1)
    assert do_match(f"{before}({constraint}){after}", text) == expected


@pytest.mark.parametrize(
    ("before", "text", "after", "expected"),
    [
        ("prefix", "prefix/", "", True),
        ("", "prefix/", "", False),
        ("", "", "", True),
        ("", "more", "more", True),
        ("prefix", "prefixmore", "more", False),
        ("prefix", "prefix/more", "more", True),
        ("prefix/", "prefix/more", "more", True),
        ("prefix", "prefix//more", "more", False),
        ("", "/more", "more", False),
        ("/", "/more", "more", True),
        ("", "/", "", False),
        ("/", "/", "", True),
        ("/", "//", "", True),
    ],
)
def test_d_wildcard_matching(before: str, after: str, text: str, expected: bool):
    """Test __d__ wildcard matching behavior."""
    _, constraint = SnakemakeWildcards.d.split(",", 1)
    assert do_match(f"{before}({constraint}){after}", text) == expected


@pytest.mark.parametrize(
    ("before", "text", "after", "expected"),
    [
        ("", "", "", True),
        ("", "more", "more", True),
        ("prefix", "prefix", "", True),
        ("prefix", "prefixmore", "more", False),
        ("prefix", "prefix", "", False),
        ("prefix", "prefix_more", "more", True),
        ("prefix", "prefix/more", "more", False),
        ("prefix/", "prefix/more", "more", True),
        ("prefix/", "prefix/_more", "more", False),
        ("prefix", "prefix/_more", "more", False),
        ("", "_more", "more", True),
        ("", "_m", "", False),
        ("", "_", "", False),
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
        ("", "", "", True),
        ("", "more", "more", True),
        ("prefix", "prefix", "", True),
        ("prefix", "prefixmore", "more", True),
        ("", "anat01/", "/", True),
        ("", "anat/", "", False),
        ("", "anat", "", False),
        ("", "anat_/", "/", False),
        ("", "anat-01/", "/", False),
        ("", "an/a/", "/", False),
        ("", "ana\n/", "/", False),
        ("", "/ana-01/", "/", False),
        ("/", "/anat01/", "", False),
        ("/", "/anat01/", "/", True),
        ("prefix", "prefixanat01/", "/", False),
        ("prefix", "prefix/anat/", "/", False),
        ("prefix/", "prefix/anat/", "", False),
        ("prefix/", "prefix/anat/", "/", True),
        ("prefix/", "prefix/anat01/more", "/more", True),
    ],
)
def test_datatype_wildcard_matching(before: str, after: str, text: str, expected: bool):
    _, constraint = SnakemakeWildcards.datatype.split(",", 1)
    assert do_match(f"{before}({constraint}){after}", text) == expected


@pytest.mark.parametrize(
    ("before", "text", "after", "expected"),
    [
        ("", "", "", True),
        ("", "more", "more", True),
        ("prefix", "prefix", "", True),
        ("prefix", "prefixmore", "more", True),
        ("", "run-0", "0", True),
        ("", "nur-0", "0", False),
        ("", "run-", "", False),
        ("", "run-0", "", False),
        ("", "run0", "0", False),
        ("", "run", "", False),
        ("", "run-_", "_", False),
        ("", "run--", "-", False),
        ("", "run-\n", "\n", False),
        ("", "run-/", "/", False),
        ("", "/run-0", "0", False),
        ("/", "/run-0", "0", True),
        ("_", "_run-0", "0", False),
        ("", "_run-0", "0", False),
        ("prefix_", "prefix_run-0", "0", False),
        ("prefix", "prefix_run-0", "0", True),
        ("prefix", "prefix/run-0", "0", False),
        ("prefix/", "prefix/run-0", "0", True),
        ("prefix/", "prefix/_run-0", "0", False),
        ("prefix/_", "prefix/_run-0", "0", False),
        ("prefix/", "prefix/run-0/more", "0/more", True),
    ],
)
def test_dummy_wildcard_matching(before: str, after: str, text: str, expected: bool):
    _, constraint = SnakemakeWildcards("run").dummy.split(",", 1)
    assert do_match(f"{before}({constraint}){after}", text) == expected


@pytest.mark.parametrize(
    ("before", "text", "after", "expected"),
    [
        ("", "", "", True),
        ("", "more", "more", True),
        ("prefix", "prefix", "", True),
        ("prefix", "prefixmore", "more", True),
        ("run-", "run-value0", "", True),
        ("", "run-value0", "", False),
        ("run-", "nur-value0", "", False),
        ("run-", "run-", "", True),
        ("run", "runvalue0", "", False),
        ("run-", "run-_", "", False),
        ("run-", "run--", "", False),
        ("run-", "run-\nfwaef", "", False),
        ("run-", "run-/", "", False),
        ("/run-", "/run-value0", "", True),
        ("prefix/run-", "prefix/run-0", "", True),
        ("prefixrun-", "prefixrun-0", "", True),
        ("run-", "run-00", "0", False),
        ("run-", "run-0more", "more", False),
        ("run-", "run-0_more", "_more", True),
        ("run-", "run-0/more", "/more", True),
        ("run-", "run-0_more", "more", False),
        ("run-", "run-0/more", "more", False),
    ],
)
def test_variable_wildcard_matching(before: str, after: str, text: str, expected: bool):
    _, constraint = SnakemakeWildcards("run").variable.split(",", 1)
    assert do_match(f"{before}({constraint}){after}", text) == expected


class TestSuffixWildcard:
    """Tests for suffix special wildcard."""

    @pytest.mark.parametrize(
        ("text", "expected"),
        [
            ("bold", "bold"),
            ("/bold", "bold"),
            ("_bold", "bold"),
            ("bol_d", "bol"),
            ("bol/d", "bol"),
        ],
    )
    def test_suffix_wildcard_matching(self, text: str, expected: str):
        """Test suffix wildcard matching behavior."""
        _, constraint = SnakemakeWildcards.suffix.split(",", 1)
        pattern = re.compile(constraint)
        result = find_non_empty_match(pattern, text)
        assert result == expected


class TestExtensionWildcard:
    """Tests for extension special wildcard."""

    @pytest.mark.parametrize(
        ("text", "expected"),
        [
            ("file.nii.gz", ".nii.gz"),
            ("filegz", ""),
            ("file.ni_i", ""),
            ("file.ni/i", ""),
            ("file.ni-i", ""),
        ],
    )
    def test_extension_wildcard_matching(self, text: str, expected: str):
        """Test extension wildcard matching behavior."""
        _, constraint = SnakemakeWildcards.extension.split(",", 1)
        pattern = re.compile(constraint)
        result = find_non_empty_match(pattern, text)
        assert result == expected


class TestSubjectSessionReplacement:
    """Test that sub and ses are replaced with subject and session."""

    @pytest.mark.parametrize(
        ("tag", "expected_name"),
        [
            ("sub", "subject"),
            ("ses", "session"),
            ("run", "run"),
        ],
    )
    def test_tag_replacement(self, tag: str, expected_name: str):
        """Test tag to wildcard name conversion."""
        wc = SnakemakeWildcards(tag)
        assert expected_name in wc.variable
        assert expected_name in wc.dummy
        assert expected_name in wc.directory
