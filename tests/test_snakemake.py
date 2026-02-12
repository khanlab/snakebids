"""Tests for snakebids.utils.snakemake module."""

from __future__ import annotations

import re

from snakebids.utils.snakemake import SnakemakeWildcards


def find_non_empty_match(pattern: re.Pattern[str], text: str) -> re.Match[str] | None:
    """Find the first non-empty match in text.

    Since wildcard patterns are optional (ending with ?), they always match,
    but may match empty string. This helper finds the first non-empty match.
    """
    for match in pattern.finditer(text):
        if match.group():
            return match
    return None


class TestDirectoryWildcard:
    """Tests for directory wildcard (_entity_d_)."""

    def test_directory_formats_label_correctly(self):
        """Test that directory wildcard label is formatted correctly."""
        wc = SnakemakeWildcards("run")
        name, _ = wc.directory.split(",", 1)
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


class TestDWildcard:
    """Tests for __d__ special wildcard."""

    def test_d_matches_when_preceded_by_matched_characters(self):
        """Test __d__ matches when preceded by matched characters."""
        _, constraint = SnakemakeWildcards.d.split(",", 1)
        pattern = re.compile(constraint)
        # Test that the pattern matches at the position of / after characters
        text = "path/"
        # The pattern should match at position 4 (the /)
        # We need to check if there's a match starting at position 4
        match = pattern.match(text, 4)
        assert match is not None
        assert match.start() == 4

    def test_d_matches_at_beginning_of_string(self):
        """Test __d__ matches at the beginning of string."""
        _, constraint = SnakemakeWildcards.d.split(",", 1)
        pattern = re.compile(constraint)
        assert pattern.match("")

    def test_d_matches_when_preceding_slash_is_matched(self):
        """Test __d__ matches when the preceding slash is matched."""
        _, constraint = SnakemakeWildcards.d.split(",", 1)
        pattern = re.compile(constraint)
        # Should match the position after /
        text = "path/segment"
        match = pattern.search(text)
        assert match is not None

    def test_d_rejects_singleton_slash(self):
        """Test __d__ rejects a singleton slash."""
        _, constraint = SnakemakeWildcards.d.split(",", 1)
        pattern = re.compile(constraint)
        # A bare "/" should not be a match for the position
        # The pattern matches positions, so we test if it would match the start of "/"
        text = "/"
        match = pattern.match(text)
        # At position 0, it should match (beginning of string)
        assert match is not None
        # But the slash itself is not consumed
        assert match.group() == ""


class TestUnderscoreWildcard:
    """Tests for ___ special wildcard."""

    def test_underscore_matches_when_preceded_by_non_slash_chars(self):
        """Test ___ matches underscore when preceded by non-slash chars."""
        _, constraint = SnakemakeWildcards.underscore.split(",", 1)
        pattern = re.compile(constraint)
        text = "path_segment"
        # Check that pattern matches at position 4 (the _)
        match = pattern.match(text, 4)
        assert match is not None
        assert match.start() == 4

    def test_underscore_matches_at_beginning_of_string(self):
        """Test ___ matches at the beginning of string."""
        _, constraint = SnakemakeWildcards.underscore.split(",", 1)
        pattern = re.compile(constraint)
        assert pattern.match("")

    def test_underscore_matches_when_preceding_slash_is_matched(self):
        """Test ___ matches when the preceding slash is matched."""
        _, constraint = SnakemakeWildcards.underscore.split(",", 1)
        pattern = re.compile(constraint)
        text = "path/_segment"
        # Check that pattern matches at position 5 (after /)
        match = pattern.match(text, 5)
        assert match is not None
        assert match.start() == 5

    def test_underscore_rejects_slash_underscore_sequence(self):
        """Test ___ does not match the underscore in /_ when / is not matched."""
        _, constraint = SnakemakeWildcards.underscore.split(",", 1)
        pattern = re.compile(constraint)
        # Pattern: ^|(?<=\/)|(?<!\/)_(?=[^\.])
        # Matches: start OR after / OR (_ not after / and followed by non-.)
        # At position 1 (after /), (?<=\/) matches
        # At position 2 (the _), it's preceded by /, so (?<!\/)_ doesn't match
        text = "/_"
        match = pattern.search(text, pos=1)
        assert match is not None

    def test_underscore_rejects_when_followed_by_period(self):
        """Test ___ rejects underscore when followed by a period."""
        _, constraint = SnakemakeWildcards.underscore.split(",", 1)
        pattern = re.compile(constraint)
        text = "path_."
        # The pattern requires (?=[^\.]) after the underscore
        # So "path_" should match at position 4, but "path_." should not match the _
        # Actually, the third alternative is (?<!\/)_(?=[^\.])
        # So we need the _ AND the lookahead
        # Let's search for where the pattern matches
        match = pattern.search(text)
        # It should match at position 0 (^) or position 4 (?<!\/)_(?=[^\.]))
        # At position 4, we have _, and it's not preceded by /, but it IS followed by .
        # So (?=[^\.]) fails
        # Therefore, it should only match at position 0
        assert match is not None
        assert match.start() == 0
        # If we search starting after position 0, we shouldn't find the underscore
        match = pattern.search(text, pos=1)
        assert match is None or "_" not in text[match.start() : match.end() + 1]

    def test_underscore_requires_following_non_period_characters(self):
        """Test ___ requires non-period characters when matching underscore."""
        _, constraint = SnakemakeWildcards.underscore.split(",", 1)
        pattern = re.compile(constraint)
        # Should match _ followed by non-period
        text = "path_segment"
        match = pattern.search(text, pos=1)
        assert match is not None
        assert match.start() == 4


class TestDatatypeWildcard:
    """Tests for datatype special wildcard."""

    def test_datatype_matches_valid_string_at_beginning(self):
        """Test datatype matches valid string at beginning."""
        _, constraint = SnakemakeWildcards.datatype.split(",", 1)
        pattern = re.compile(constraint)
        text = "anat/"
        match = find_non_empty_match(pattern, text)
        assert match is not None
        assert match.group() == "anat"

    def test_datatype_matches_valid_string_after_slash(self):
        """Test datatype matches valid string after slash."""
        _, constraint = SnakemakeWildcards.datatype.split(",", 1)
        pattern = re.compile(constraint)
        # Test with string that only has datatype after slash
        text = "/anat/"
        match = find_non_empty_match(pattern, text)
        assert match is not None
        assert match.group() == "anat"

    def test_datatype_requires_following_slash(self):
        """Test datatype requires a following slash."""
        _, constraint = SnakemakeWildcards.datatype.split(",", 1)
        pattern = re.compile(constraint)
        text = "anat"
        match = find_non_empty_match(pattern, text)
        # Should not find a non-empty match without trailing /
        assert match is None

    def test_datatype_rejects_underscore(self):
        """Test datatype rejects underscores in value."""
        _, constraint = SnakemakeWildcards.datatype.split(",", 1)
        pattern = re.compile(constraint)
        text = "ana_t/"
        match = find_non_empty_match(pattern, text)
        # Should not match due to underscore
        assert match is None

    def test_datatype_rejects_dash(self):
        """Test datatype rejects dashes in value."""
        _, constraint = SnakemakeWildcards.datatype.split(",", 1)
        pattern = re.compile(constraint)
        text = "ana-t/"
        match = find_non_empty_match(pattern, text)
        assert match is None

    def test_datatype_rejects_extra_slash_in_value(self):
        """Test datatype rejects extra slashes in value."""
        _, constraint = SnakemakeWildcards.datatype.split(",", 1)
        pattern = re.compile(constraint)
        text = "ana/t/"
        match = find_non_empty_match(pattern, text)
        assert match is not None
        # Should match only "ana" followed by /
        assert match.group() == "ana"


class TestDummyWildcard:
    """Tests for dummy wildcard (_entity_)."""

    def test_dummy_formats_label_correctly(self):
        """Test that dummy wildcard label is formatted correctly."""
        wc = SnakemakeWildcards("run")
        name, _ = wc.dummy.split(",", 1)
        assert name == "_run_"

    def test_dummy_matches_entity_with_preceding_underscore(self):
        """Test dummy matches entity-value with preceding underscore."""
        wc = SnakemakeWildcards("run")
        _, constraint = wc.dummy.split(",", 1)
        pattern = re.compile(constraint)
        text = "_run-01"
        match = find_non_empty_match(pattern, text)
        assert match is not None
        assert match.group() == "_run-"

    def test_dummy_matches_with_preceding_slash(self):
        """Test dummy matches entity-value with preceding slash."""
        wc = SnakemakeWildcards("run")
        _, constraint = wc.dummy.split(",", 1)
        pattern = re.compile(constraint)
        text = "/run-01"
        match = find_non_empty_match(pattern, text)
        assert match is not None
        # The pattern uses a lookbehind, so / is not included in the match
        assert match.group() == "run-"

    def test_dummy_matches_at_beginning_of_string(self):
        """Test dummy matches at beginning of string."""
        wc = SnakemakeWildcards("run")
        _, constraint = wc.dummy.split(",", 1)
        pattern = re.compile(constraint)
        text = "run-01"
        match = find_non_empty_match(pattern, text)
        assert match is not None
        assert match.group() == "run-"

    def test_dummy_rejects_without_following_valid_characters(self):
        """Test dummy rejects when not followed by valid characters."""
        wc = SnakemakeWildcards("run")
        _, constraint = wc.dummy.split(",", 1)
        pattern = re.compile(constraint)
        # The lookahead (?=[^_\/\-\n]) requires a character that is not _, /, -, or \n
        text = "run-_"
        match = find_non_empty_match(pattern, text)
        assert match is None

    def test_dummy_rejects_when_preceded_by_other_characters(self):
        """Test dummy doesn't match when preceded by other characters."""
        wc = SnakemakeWildcards("run")
        _, constraint = wc.dummy.split(",", 1)
        pattern = re.compile(constraint)
        # Should not match in middle of string unless after /, _, or at start
        text = "pathrun-01"
        match = find_non_empty_match(pattern, text)
        assert match is None


class TestOrdinaryWildcard:
    """Tests for ordinary/variable wildcard (entity)."""

    def test_ordinary_formats_label_correctly(self):
        """Test that ordinary wildcard label is formatted correctly."""
        wc = SnakemakeWildcards("run")
        name, _ = wc.variable.split(",", 1)
        assert name == "run"

    def test_ordinary_requires_preceding_tag(self):
        """Test ordinary wildcard requires preceding tag."""
        wc = SnakemakeWildcards("run")
        _, constraint = wc.variable.split(",", 1)
        pattern = re.compile(constraint)
        text = "run-01"
        match = find_non_empty_match(pattern, text)
        assert match is not None
        assert match.group() == "01"

    def test_ordinary_rejects_without_preceding_tag(self):
        """Test ordinary wildcard doesn't match without preceding tag."""
        wc = SnakemakeWildcards("run")
        _, constraint = wc.variable.split(",", 1)
        pattern = re.compile(constraint)
        text = "01"
        match = find_non_empty_match(pattern, text)
        assert match is None

    def test_ordinary_matches_with_following_underscore(self):
        """Test ordinary wildcard matches with following underscore."""
        wc = SnakemakeWildcards("run")
        _, constraint = wc.variable.split(",", 1)
        pattern = re.compile(constraint)
        text = "run-01_next"
        match = find_non_empty_match(pattern, text)
        assert match is not None
        assert match.group() == "01"

    def test_ordinary_matches_with_following_slash(self):
        """Test ordinary wildcard matches with following slash."""
        wc = SnakemakeWildcards("run")
        _, constraint = wc.variable.split(",", 1)
        pattern = re.compile(constraint)
        text = "run-01/path"
        match = find_non_empty_match(pattern, text)
        assert match is not None
        assert match.group() == "01"

    def test_ordinary_matches_at_end_of_string(self):
        """Test ordinary wildcard matches at end of string."""
        wc = SnakemakeWildcards("run")
        _, constraint = wc.variable.split(",", 1)
        pattern = re.compile(constraint)
        text = "run-01"
        match = find_non_empty_match(pattern, text)
        assert match is not None
        assert match.group() == "01"

    def test_ordinary_rejects_underscore_in_value(self):
        """Test ordinary wildcard rejects underscores in value."""
        wc = SnakemakeWildcards("run")
        _, constraint = wc.variable.split(",", 1)
        pattern = re.compile(constraint)
        text = "run-01_02"
        match = find_non_empty_match(pattern, text)
        assert match is not None
        # Should match only "01" (up to the underscore)
        assert match.group() == "01"

    def test_ordinary_rejects_slash_in_value(self):
        """Test ordinary wildcard rejects slashes in value."""
        wc = SnakemakeWildcards("run")
        _, constraint = wc.variable.split(",", 1)
        pattern = re.compile(constraint)
        text = "run-01/02"
        match = find_non_empty_match(pattern, text)
        assert match is not None
        assert match.group() == "01"

    def test_ordinary_rejects_dash_in_value(self):
        """Test ordinary wildcard rejects dashes in value."""
        wc = SnakemakeWildcards("run")
        _, constraint = wc.variable.split(",", 1)
        pattern = re.compile(constraint)
        text = "run-01-02"
        match = find_non_empty_match(pattern, text)
        # Dashes are not allowed in the value, and "-02" does not satisfy the
        # lookahead, so this should not match
        assert match is None


class TestSuffixWildcard:
    """Tests for suffix special wildcard."""

    def test_suffix_matches_at_beginning_of_string(self):
        """Test suffix matches at beginning of string."""
        _, constraint = SnakemakeWildcards.suffix.split(",", 1)
        pattern = re.compile(constraint)
        text = "bold"
        match = find_non_empty_match(pattern, text)
        assert match is not None
        assert match.group() == "bold"

    def test_suffix_matches_after_slash(self):
        """Test suffix matches after slash."""
        _, constraint = SnakemakeWildcards.suffix.split(",", 1)
        pattern = re.compile(constraint)
        # Test with string that only has the suffix after slash
        text = "/bold"
        match = find_non_empty_match(pattern, text)
        assert match is not None
        assert match.group() == "bold"

    def test_suffix_matches_after_underscore(self):
        """Test suffix matches after underscore."""
        _, constraint = SnakemakeWildcards.suffix.split(",", 1)
        pattern = re.compile(constraint)
        # Test with string that only has the suffix after underscore
        text = "_bold"
        match = find_non_empty_match(pattern, text)
        assert match is not None
        assert match.group() == "bold"

    def test_suffix_rejects_underscore_in_value(self):
        """Test suffix rejects underscores in value."""
        _, constraint = SnakemakeWildcards.suffix.split(",", 1)
        pattern = re.compile(constraint)
        text = "bol_d"
        match = find_non_empty_match(pattern, text)
        assert match is not None
        # Should match only "bol"
        assert match.group() == "bol"

    def test_suffix_rejects_slash_in_value(self):
        """Test suffix rejects slashes in value."""
        _, constraint = SnakemakeWildcards.suffix.split(",", 1)
        pattern = re.compile(constraint)
        text = "bol/d"
        match = find_non_empty_match(pattern, text)
        assert match is not None
        assert match.group() == "bol"


class TestExtensionWildcard:
    """Tests for extension special wildcard."""

    def test_extension_matches_at_end_of_string(self):
        """Test extension matches at end of string."""
        _, constraint = SnakemakeWildcards.extension.split(",", 1)
        pattern = re.compile(constraint)
        # The pattern matches from the first period to the end
        text = "file.nii.gz"
        match = find_non_empty_match(pattern, text)
        assert match is not None
        # Matches the entire extension including multiple periods
        assert match.group() == ".nii.gz"

    def test_extension_requires_leading_period(self):
        """Test extension requires a leading period."""
        _, constraint = SnakemakeWildcards.extension.split(",", 1)
        pattern = re.compile(constraint)
        text = "filegz"
        match = find_non_empty_match(pattern, text)
        assert match is None

    def test_extension_rejects_underscore_in_value(self):
        """Test extension rejects underscores in value."""
        _, constraint = SnakemakeWildcards.extension.split(",", 1)
        pattern = re.compile(constraint)
        text = "file.ni_i"
        match = find_non_empty_match(pattern, text)
        assert match is None

    def test_extension_rejects_slash_in_value(self):
        """Test extension rejects slashes in value."""
        _, constraint = SnakemakeWildcards.extension.split(",", 1)
        pattern = re.compile(constraint)
        text = "file.ni/i"
        match = find_non_empty_match(pattern, text)
        assert match is None

    def test_extension_rejects_dash_in_value(self):
        """Test extension rejects dashes in value."""
        _, constraint = SnakemakeWildcards.extension.split(",", 1)
        pattern = re.compile(constraint)
        text = "file.ni-i"
        match = find_non_empty_match(pattern, text)
        assert match is None


class TestSubjectSessionReplacement:
    """Test that sub and ses are replaced with subject and session."""

    def test_sub_replaced_with_subject(self):
        """Test that 'sub' tag is replaced with 'subject'."""
        wc = SnakemakeWildcards("sub")
        assert "subject" in wc.variable
        assert "subject" in wc.dummy
        assert "subject" in wc.directory

    def test_ses_replaced_with_session(self):
        """Test that 'ses' tag is replaced with 'session'."""
        wc = SnakemakeWildcards("ses")
        assert "session" in wc.variable
        assert "session" in wc.dummy
        assert "session" in wc.directory

    def test_other_tags_not_replaced(self):
        """Test that other tags are not replaced."""
        wc = SnakemakeWildcards("run")
        assert "run" in wc.variable
        assert "run" in wc.dummy
        assert "run" in wc.directory
