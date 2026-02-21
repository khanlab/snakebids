"""Tests for snakebids.utils.snakemake_templates module."""

from __future__ import annotations

import re

import pytest
from hypothesis import given
from hypothesis import strategies as st

from snakebids.utils.snakemake_templates import SnakemakeFormatter, SnakemakeWildcards
from tests.test_snakemake_templates.strategies import safe_field_names


def do_match(pattern: str, text: str) -> str | None:
    match = re.match(f"(?:{pattern})$", text)
    if match is not None:
        return match.group(1)
    return None


class TestCorrectLabels:
    @given(entity=safe_field_names())
    def test_directory_formats_label_correctly(self, entity: str):
        """Test that directory wildcard label is formatted correctly."""
        wc = SnakemakeWildcards(entity)
        formatter = SnakemakeFormatter()
        assert (
            formatter.format(str(wc.directory), **{wc._wildcard: "foo"})
            == f"{wc._tag}-foo/"
        )

    @given(entity=safe_field_names())
    def test_dummy_formats_label_correctly(self, entity: str):
        """Test that directory wildcard label is formatted correctly."""
        wc = SnakemakeWildcards(entity)
        formatter = SnakemakeFormatter()
        assert formatter.format(str(wc.dummy), **{wc._wildcard: "foo"}) == f"{wc._tag}-"

    @given(entity=safe_field_names(min_size=1).filter(lambda s: not s.isdigit()))
    def test_variable_formats_label_correctly(self, entity: str):
        """Test that directory wildcard label is formatted correctly."""
        wc = SnakemakeWildcards(entity)
        formatter = SnakemakeFormatter()
        assert formatter.format(str(wc.variable), **{wc._wildcard: "foo"}) == "foo"

    def test_underscore_formats_label_correctly(self):
        """Test that directory wildcard label is formatted correctly."""
        formatter = SnakemakeFormatter()
        assert formatter.format(f"foo{SnakemakeWildcards.underscore}") == "foo_"

    def test_slash_formats_label_correctly(self):
        """Test that directory wildcard label is formatted correctly."""
        formatter = SnakemakeFormatter()
        assert (
            formatter.format(f"anat{SnakemakeWildcards.slash}", datatype="anat")
            == "anat/"
        )


@pytest.mark.parametrize("entity", ["datatype", "suffix", "extension"])
def test_special_entities_formats_label_correctly(entity: str):
    """Test that directory wildcard label is formatted correctly."""
    formatter = SnakemakeFormatter()
    assert (
        formatter.format(f"{getattr(SnakemakeWildcards, entity)}", **{entity: "foo"})
        == "foo"
    )


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
        ("", "run-\n/", "", "run-\n/"),
        ("", "run-./", "", None),
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
def test_directory_wildcard_matching(
    before: str, text: str, after: str, expected: str | None
):
    """Test directory wildcard matching behavior."""
    wc = SnakemakeWildcards("run")
    constraint = wc.directory.constraint
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
def test_d_wildcard_matching(before: str, text: str, after: str, expected: str | None):
    """Test __d__ wildcard matching behavior."""
    constraint = SnakemakeWildcards.slash.constraint
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
    before: str, text: str, after: str, expected: str | None
):
    constraint = SnakemakeWildcards.underscore.constraint
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
        ("", "ana\n/", "/", "ana\n"),
        ("", "anat./", "/", None),
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
def test_datatype_wildcard_matching(
    before: str, text: str, after: str, expected: str | None
):
    constraint = SnakemakeWildcards.datatype.constraint
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
        ("", "run-\n", "\n", "run-"),
        ("", "run-.", ".", None),
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
def test_dummy_wildcard_matching(
    before: str, text: str, after: str, expected: str | None
):
    constraint = SnakemakeWildcards("run").dummy.constraint
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
        ("run-", "run-\nfwaef", "", "\nfwaef"),
        ("run-", "run-.value", "", None),
        ("run-", "run-value.ext", "", None),
        ("run-", "run-value.ext", "ext", None),
        ("run-", "run-value.ext", ".ext", "value"),
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
def test_ordinary_wildcard_matching(
    before: str, text: str, after: str, expected: str | None
):
    constraint = SnakemakeWildcards("run").variable.constraint
    assert do_match(f"{before}({constraint}){after}", text) == expected


@pytest.mark.parametrize(
    ("before", "text", "after", "expected"),
    [
        ("", "", "", ""),
        ("", "more", "more", ""),
        ("prefix", "prefix", "", ""),
        ("prefix", "prefixmore", "more", ""),
        ("", "suffix.ext", "", None),
        ("", "suffix-ext", "", None),
        ("", "suffix_ext", "", None),
        ("", "suffix/ext", "", None),
        ("", "_suffix", "", None),
        ("_", "_suffix", "", "suffix"),
        ("/_", "/_suffix", "", "suffix"),
        ("", "/suffix", "", None),
        ("/", "/suffix", "", "suffix"),
        ("/", "/suffixmore", "more", "suffix"),
        ("/", "/suffix_more", "_more", "suffix"),
        ("/", "/suffix/more", "/more", "suffix"),
        ("/", "/suffix-more", "-more", "suffix"),
        ("/", "/suffix.more", ".more", "suffix"),
    ],
)
def test_suffix_wildcard_matching(
    before: str, text: str, after: str, expected: str | None
):
    constraint = SnakemakeWildcards.suffix.constraint
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
    before: str, text: str, after: str, expected: str | None
):
    constraint = SnakemakeWildcards.extension.constraint
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
    assert expected_name in wc.variable.label
    assert expected_name in wc.dummy.label
    assert expected_name in wc.directory.label


@pytest.mark.parametrize(
    ("tag", "wildcard_label"),
    [
        ("sub", "subject"),
        ("ses", "session"),
    ],
)
def test_subject_session_label_vs_tag_in_constraints(tag: str, wildcard_label: str):
    """Regression test: subject/session wildcard labels and constraint tag prefixes.

    The wildcard *label* (used inside ``{...}``) must remain the full name
    ("subject"/"session"), while the regex *constraints* must reference the
    short BIDS tag prefix ("sub-"/"ses-") that actually appears in file paths.
    These two values must not be confused.
    """
    wc = SnakemakeWildcards(tag)

    # Wildcard labels must use the full wildcard name, and constraints must
    # reference the short BIDS tag prefix, not the full wildcard label
    for wildcard in (wc.variable, wc.dummy, wc.directory):
        # Wildcard labels must use the full wildcard name
        assert wildcard_label in wildcard.label

        constraint = wildcard.constraint
        # The regex constraints use escaped hyphens (e.g. "sub\-")
        assert f"{tag}\\-" in constraint
        # Ensure constraints do NOT use the full wildcard label as tag prefix
        assert f"{wildcard_label}\\-" not in constraint


class TestBraceWrappedOutput:
    @pytest.mark.parametrize("attribute", ["variable", "dummy", "directory"])
    @given(entity=st.text())
    def test_dynamic_wildcards_are_brace_wrapped(self, entity: str, attribute: str):
        """Test that variable wildcard returns a brace-wrapped string."""
        wildcard = getattr(SnakemakeWildcards(entity), attribute).wildcard
        assert wildcard.startswith("{")
        assert wildcard.endswith("}")

    @pytest.mark.parametrize(
        "attr",
        ["underscore", "slash", "datatype", "suffix", "extension"],
    )
    def test_class_attributes_are_brace_wrapped(self, attr: str):
        value = getattr(SnakemakeWildcards, attr).wildcard
        assert value.startswith("{")
        assert value.endswith("}")


class TestConstraintsAccessor:
    """Test that .constraint exposes raw regex constraint strings."""

    @given(entity=st.text())
    def test_wildcard_wraps_constraint(self, entity: str):
        wc = SnakemakeWildcards(entity)
        assert wc.variable.constraint in wc.variable.wildcard
        assert wc.dummy.constraint in wc.dummy.wildcard
        assert wc.directory.constraint in wc.directory.wildcard

    @pytest.mark.parametrize(
        "attribute",
        [
            "underscore",
            "slash",
            "datatype",
            "suffix",
            "extension",
        ],
    )
    def test_class_wildcard_wraps_constraint(self, attribute: str):
        """Test that class-level brace-wrapped attr contains the constraint."""
        wildcard = getattr(SnakemakeWildcards, attribute)
        assert wildcard.constraint in wildcard.wildcard
