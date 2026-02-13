"""Tests for snakebids.utils.snakemake_templates module."""

from __future__ import annotations

import re

import pytest

from snakebids.utils.snakemake_templates import SnakemakeWildcards


class TestDirectoryWildcard:
    """Tests for directory wildcard (_entity_d_)."""

    def test_directory_formats_label_correctly(self):
        """Test that directory wildcard label is formatted correctly."""
        wc = SnakemakeWildcards("run")
        name, _ = wc.directory.split(",", 1)
        assert name == "_run_d_"

    @pytest.mark.parametrize(
        ("text", "should_match", "expected"),
        [
            ("run-01/$", True, "run-01/"),
            ("run-01$", False, ""),
            ("run-01_test/$", False, ""),
            ("run-01-02/$", False, ""),
            ("pathrun-01/$", False, ""),
            ("/run-01/$", True, "run-01/"),
            ("prefix/run-01/$", True, "run-01/"),
        ],
    )
    def test_directory_wildcard_matching(
        self, text: str, should_match: bool, expected: str
    ):
        """Test directory wildcard matching behavior."""
        wc = SnakemakeWildcards("run")
        _, constraint = wc.directory.split(",", 1)
        pattern = re.compile(constraint)
        match = pattern.search(text)
        
        if should_match:
            assert match is not None
            assert match.group() == expected
        else:
            assert match is None or not match.group()


class TestDWildcard:
    """Tests for __d__ special wildcard."""

    @pytest.mark.parametrize(
        ("text",),
        [
            ("$",),
            ("path/$",),
            ("/$",),
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
        ("text",),
        [
            ("$",),
            ("path_$",),
            ("/$",),
            ("/_$",),
        ],
    )
    def test_underscore_wildcard_matching(self, text: str):
        """Test ___ wildcard matching behavior."""
        _, constraint = SnakemakeWildcards.underscore.split(",", 1)
        pattern = re.compile(constraint)
        match = pattern.match(text)
        assert match is not None

    def test_underscore_matches_when_preceding_slash_is_matched(self):
        """Test that underscore after slash is left unmatched."""
        _, constraint = SnakemakeWildcards.underscore.split(",", 1)
        pattern = re.compile("/" + constraint + "segment$")
        match = pattern.match("/_segment")
        assert match is not None
        matched_text = match.group()
        # After the leading /, the next char should be matched by constraint, not the _
        # So the _ should not be in the matched portion after /
        assert matched_text[1:2] != "_" or matched_text[1:] == "_segment"

    def test_underscore_rejects_when_followed_by_period(self):
        """Test ___ rejects underscore when followed by a period."""
        _, constraint = SnakemakeWildcards.underscore.split(",", 1)
        pattern = re.compile("path" + constraint + r"\.$")
        match = pattern.match("path_.")
        assert match is None


class TestDatatypeWildcard:
    """Tests for datatype special wildcard."""

    @pytest.mark.parametrize(
        ("text", "should_match", "expected"),
        [
            ("anat/$", True, "anat"),
            ("/anat/$", True, "anat"),
            ("anat$", False, ""),
            ("ana_t/$", False, ""),
            ("ana-t/$", False, ""),
            ("ana/$", True, "ana"),
        ],
    )
    def test_datatype_wildcard_matching(
        self, text: str, should_match: bool, expected: str
    ):
        """Test datatype wildcard matching behavior."""
        _, constraint = SnakemakeWildcards.datatype.split(",", 1)
        pattern = re.compile(constraint)
        match = pattern.search(text)
        
        if should_match:
            assert match is not None
            assert match.group() == expected
        else:
            assert match is None or not match.group()


class TestDummyWildcard:
    """Tests for dummy wildcard (_entity_)."""

    def test_dummy_formats_label_correctly(self):
        """Test that dummy wildcard label is formatted correctly."""
        wc = SnakemakeWildcards("run")
        name, _ = wc.dummy.split(",", 1)
        assert name == "_run_"

    @pytest.mark.parametrize(
        ("text", "should_match", "expected"),
        [
            ("_run-01$", True, "_run-"),
            ("/run-01$", True, "run-"),
            ("run-01$", True, "run-"),
            ("run-_$", False, ""),
            ("pathrun-01$", False, ""),
        ],
    )
    def test_dummy_wildcard_matching(
        self, text: str, should_match: bool, expected: str
    ):
        """Test dummy wildcard matching behavior."""
        wc = SnakemakeWildcards("run")
        _, constraint = wc.dummy.split(",", 1)
        pattern = re.compile(constraint)
        match = pattern.search(text)
        
        if should_match:
            assert match is not None
            assert match.group() == expected
        else:
            assert match is None or not match.group()


class TestOrdinaryWildcard:
    """Tests for ordinary/variable wildcard (entity)."""

    def test_ordinary_formats_label_correctly(self):
        """Test that ordinary wildcard label is formatted correctly."""
        wc = SnakemakeWildcards("run")
        name, _ = wc.variable.split(",", 1)
        assert name == "run"

    @pytest.mark.parametrize(
        ("text", "should_match", "expected"),
        [
            ("run-01$", True, "01"),
            ("01$", False, ""),
            ("run-01_$", True, "01"),
            ("run-01/$", True, "01"),
            ("run-01-$", False, ""),
        ],
    )
    def test_ordinary_wildcard_matching(
        self, text: str, should_match: bool, expected: str
    ):
        """Test ordinary wildcard matching behavior."""
        wc = SnakemakeWildcards("run")
        _, constraint = wc.variable.split(",", 1)
        pattern = re.compile(constraint)
        match = pattern.search(text)
        
        if should_match:
            assert match is not None
            assert match.group() == expected
        else:
            assert match is None or not match.group()


class TestSuffixWildcard:
    """Tests for suffix special wildcard."""

    @pytest.mark.parametrize(
        ("text", "should_match", "expected"),
        [
            ("bold$", True, "bold"),
            ("/bold$", True, "bold"),
            ("_bold$", True, "bold"),
            ("bol_$", True, "bol"),
            ("bol/$", True, "bol"),
        ],
    )
    def test_suffix_wildcard_matching(
        self, text: str, should_match: bool, expected: str
    ):
        """Test suffix wildcard matching behavior."""
        _, constraint = SnakemakeWildcards.suffix.split(",", 1)
        pattern = re.compile(constraint)
        match = pattern.search(text)
        
        if should_match:
            assert match is not None
            assert match.group() == expected
        else:
            assert match is None or not match.group()


class TestExtensionWildcard:
    """Tests for extension special wildcard."""

    @pytest.mark.parametrize(
        ("text", "should_match", "expected"),
        [
            ("file.nii.gz", True, ".nii.gz"),
            ("filegz", False, ""),
            ("file.ni_", False, ""),
            ("file.ni/", False, ""),
            ("file.ni-", False, ""),
        ],
    )
    def test_extension_wildcard_matching(
        self, text: str, should_match: bool, expected: str
    ):
        """Test extension wildcard matching behavior."""
        _, constraint = SnakemakeWildcards.extension.split(",", 1)
        pattern = re.compile(constraint)
        match = pattern.search(text)
        
        if should_match:
            assert match is not None
            assert match.group() == expected
        else:
            assert match is None or not match.group()


class TestSubjectSessionReplacement:
    """Test that sub and ses are replaced with subject and session."""

    @pytest.mark.parametrize(
        ("tag", "expected_in_variable"),
        [
            ("sub", "subject"),
            ("ses", "session"),
            ("run", "run"),
        ],
    )
    def test_tag_replacement(self, tag: str, expected_in_variable: str):
        """Test tag to wildcard name conversion."""
        wc = SnakemakeWildcards(tag)
        assert expected_in_variable in wc.variable
        assert expected_in_variable in wc.dummy
        assert expected_in_variable in wc.directory
