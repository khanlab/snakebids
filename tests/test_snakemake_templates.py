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


class TestDirectoryWildcard:
    """Tests for directory wildcard (_entity_d_)."""

    def test_directory_formats_label_correctly(self):
        """Test that directory wildcard label is formatted correctly."""
        wc = SnakemakeWildcards("run")
        name, _ = wc.directory.split(",", 1)
        assert name == "_run_d_"

    @pytest.mark.parametrize(
        ("text", "expected"),
        [
            ("run-01/", "run-01/"),
            ("run-01", ""),
            ("run-01_test/", ""),
            ("run-01-02/", ""),
            ("pathrun-01/", ""),
            ("/run-01/rest", "run-01/"),
            ("prefix/run-01/more", "run-01/"),
            ("prefixrun-01/", ""),  # Leading chars without slash should not match
        ],
    )
    def test_directory_wildcard_matching(self, text: str, expected: str):
        """Test directory wildcard matching behavior."""
        wc = SnakemakeWildcards("run")
        _, constraint = wc.directory.split(",", 1)
        pattern = re.compile(constraint)
        result = find_non_empty_match(pattern, text)
        assert result == expected


class TestDWildcard:
    """Tests for __d__ special wildcard."""

    @pytest.mark.parametrize(
        "text",
        [
            "",
            "path/",
            "/",
        ],
    )
    def test_d_wildcard_matching(self, text: str):
        """Test __d__ wildcard matching behavior."""
        _, constraint = SnakemakeWildcards.d.split(",", 1)
        pattern = re.compile(constraint)
        match = pattern.match(text)
        assert match is not None


class TestUnderscoreWildcard:
    """Tests for ___ special wildcard."""

    @pytest.mark.parametrize(
        "text",
        [
            "",
            "path_",
            "/",
            "/_",
        ],
    )
    def test_underscore_wildcard_matching(self, text: str):
        """Test ___ wildcard matching behavior."""
        _, constraint = SnakemakeWildcards.underscore.split(",", 1)
        pattern = re.compile(constraint)
        match = pattern.match(text)
        assert match is not None

    def test_underscore_matches_when_preceding_slash_is_matched(self):
        """Test that underscore after slash doesn't include the underscore in body."""
        _, constraint = SnakemakeWildcards.underscore.split(",", 1)
        pattern = re.compile("/" + "(" + constraint + ")" + "seg")
        match = pattern.match("/seg")
        assert match is not None

    def test_underscore_rejects_when_followed_by_period(self):
        """Test ___ rejects underscore when followed by a period."""
        _, constraint = SnakemakeWildcards.underscore.split(",", 1)
        pattern = re.compile("path" + constraint + r"\.")
        match = pattern.match("path_.")
        assert match is None


class TestDatatypeWildcard:
    """Tests for datatype special wildcard."""

    @pytest.mark.parametrize(
        ("text", "expected"),
        [
            ("anat/", "anat"),
            ("/anat/more", "anat"),
            ("anat", ""),
            ("ana_t/", ""),
            ("ana-t/", ""),
            ("ana/more", "ana"),
        ],
    )
    def test_datatype_wildcard_matching(self, text: str, expected: str):
        """Test datatype wildcard matching behavior."""
        _, constraint = SnakemakeWildcards.datatype.split(",", 1)
        pattern = re.compile(constraint)
        result = find_non_empty_match(pattern, text)
        assert result == expected


class TestDummyWildcard:
    """Tests for dummy wildcard (_entity_)."""

    def test_dummy_formats_label_correctly(self):
        """Test that dummy wildcard label is formatted correctly."""
        wc = SnakemakeWildcards("run")
        name, _ = wc.dummy.split(",", 1)
        assert name == "_run_"

    @pytest.mark.parametrize(
        ("text", "expected"),
        [
            ("_run-01", "_run-"),
            ("/run-01", "run-"),
            ("run-01", "run-"),
            ("run-_", ""),
            ("pathrun-01", ""),
        ],
    )
    def test_dummy_wildcard_matching(self, text: str, expected: str):
        """Test dummy wildcard matching behavior."""
        wc = SnakemakeWildcards("run")
        _, constraint = wc.dummy.split(",", 1)
        pattern = re.compile(constraint)
        result = find_non_empty_match(pattern, text)
        assert result == expected


class TestOrdinaryWildcard:
    """Tests for ordinary/variable wildcard (entity)."""

    def test_ordinary_formats_label_correctly(self):
        """Test that ordinary wildcard label is formatted correctly."""
        wc = SnakemakeWildcards("run")
        name, _ = wc.variable.split(",", 1)
        assert name == "run"

    @pytest.mark.parametrize(
        ("text", "expected"),
        [
            ("run-01", "01"),
            ("01", ""),
            ("run-01_next", "01"),
            ("run-01/path", "01"),
            ("run-01-02", ""),
        ],
    )
    def test_ordinary_wildcard_matching(self, text: str, expected: str):
        """Test ordinary wildcard matching behavior."""
        wc = SnakemakeWildcards("run")
        _, constraint = wc.variable.split(",", 1)
        pattern = re.compile(constraint)
        result = find_non_empty_match(pattern, text)
        assert result == expected


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
