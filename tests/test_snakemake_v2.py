"""Tests for snakebids.utils.snakemake module - focusing on what patterns capture."""

from __future__ import annotations

import re

from snakebids.utils.snakemake import SnakemakeWildcards


def find_non_empty_match(pattern, text):
    """Find the first non-empty match in text."""
    for match in pattern.finditer(text):
        if match.group():
            return match
    return None


class TestDirectoryWildcard:
    """Tests for directory wildcard (_entity_d_)."""

    def test_directory_formats_label_correctly(self):
        """Test that directory wildcard label is formatted correctly."""
        wc = SnakemakeWildcards("run")
        name, constraint = wc.directory.split(",", 1)
        assert name == "_run_d_"

    def test_directory_matches_one_directory(self):
        """Test that directory wildcard matches a single directory with tag."""
        wc = SnakemakeWildcards("run")
        _, constraint = wc.directory.split(",", 1)
        pattern = re.compile(constraint)
        match = find_non_empty_match(pattern, "run-01/")
        assert match is not None
        assert match.group() == "run-01/"

    def test_directory_requires_trailing_slash(self):
        """Test that directory wildcard requires a trailing slash."""
        wc = SnakemakeWildcards("run")
        _, constraint = wc.directory.split(",", 1)
        pattern = re.compile(constraint)
        match = find_non_empty_match(pattern, "run-01")
        assert match is None

    def test_directory_rejects_underscore_in_directory(self):
        """Test that directory wildcard rejects underscores in the value."""
        wc = SnakemakeWildcards("run")
        _, constraint = wc.directory.split(",", 1)
        pattern = re.compile(constraint)
        match = find_non_empty_match(pattern, "run-01_test/")
        assert match is None

    def test_directory_rejects_multiple_dashes_in_directory(self):
        """Test that directory wildcard rejects multiple dashes."""
        wc = SnakemakeWildcards("run")
        _, constraint = wc.directory.split(",", 1)
        pattern = re.compile(constraint)
        match = find_non_empty_match(pattern, "run-01-02/")
        assert match is None

    def test_directory_rejects_umatched_leading_chars(self):
        """Test that directory wildcard rejects unmatched leading characters."""
        wc = SnakemakeWildcards("run")
        _, constraint = wc.directory.split(",", 1)
        pattern = re.compile(constraint)
        # Should not match in the middle of a path without preceding /
        match = find_non_empty_match(pattern, "pathrun-01/")
        assert match is None
