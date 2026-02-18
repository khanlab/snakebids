"""Focused unit tests for the querying module.

These tests provide faster, more targeted validation of the query module
functionality without relying on hypothesis property-based testing.
"""

from __future__ import annotations

import re

import pytest
from bids.layout import Query

from snakebids.core._querying import (
    PostFilter,
    UnifiedFilter,
    _compile_filters,
    _InvalidKeyError,
    _TooFewKeysError,
    _TooManyKeysError,
)


class TestPostFilterUnit:
    """Unit tests for PostFilter class."""

    def test_postfilter_initialized_empty(self):
        pf = PostFilter()
        assert pf.inclusions == {}
        assert pf.exclusions == {}

    @pytest.mark.parametrize(
        "other",
        [
            None,
            42,
            "string",
            [],
            {},
            object(),
        ],
    )
    def test_postfilter_not_equal_to_other_types(self, other: object):
        assert PostFilter() != other

    def test_postfilter_equality(self):
        pf1 = PostFilter()
        pf2 = PostFilter()
        pf1.add_filter("subject", "01", None)
        pf1.add_filter("session", None, "excluded")
        pf2.add_filter("subject", "01", None)
        pf2.add_filter("session", None, "excluded")
        assert pf1 == pf2

    def test_postfilter_inequality(self):
        pf1 = PostFilter()
        pf2 = PostFilter()
        pf1.add_filter("subject", "01", None)
        pf2.add_filter("subject", "02", None)
        assert pf1 != pf2

    @pytest.mark.parametrize(
        ("key", "val"),
        [
            ("subject", "01"),
            ("session", "1"),
            ("run", ["1", "2", "3"]),
            ("task", ["rest"]),
        ],
    )
    def test_add_filter_turns_inclusions_into_list(self, key: str, val: str | list[str]):
        pf = PostFilter()
        pf.add_filter(key, val, None)
        assert isinstance(pf.inclusions[key], list)
        if isinstance(val, str):
            assert pf.inclusions[key] == [val]
        else:
            assert pf.inclusions[key] == val

    @pytest.mark.parametrize(
        ("key", "vals"),
        [
            ("subject", "01"),
            ("session", "1"),
            ("run", ["1", "2"]),
            ("task", ["rest", "motor"]),
        ],
    )
    def test_add_filter_turns_exclusions_into_regex(self, key: str, vals: str | list[str]):
        pf = PostFilter()
        pf.add_filter(key, None, vals)
        assert key in pf.exclusions
        excl = pf.exclusions[key]
        assert isinstance(excl, list)
        assert len(excl) == 1
        assert isinstance(excl[0], str)

    def test_empty_list_of_exclusions_treated_as_none(self):
        pf = PostFilter()
        pf.add_filter("subject", None, iter([]))
        assert "subject" not in pf.exclusions

    @pytest.mark.parametrize(
        ("exclusion", "should_exclude", "should_include"),
        [
            ("01", ["01"], ["02", "03", "001"]),
            (["01", "02"], ["01", "02"], ["03", "001", "012"]),
            (["sub01", "sub02"], ["sub01", "sub02"], ["sub03", "sub001"]),
        ],
    )
    def test_exclusion_regex_excludes_correct_values(
        self,
        exclusion: str | list[str],
        should_exclude: list[str],
        should_include: list[str],
    ):
        pf = PostFilter()
        pf.add_filter("subject", None, exclusion)
        pattern = pf.exclusions["subject"][0]

        # Should not match excluded values
        for ex in should_exclude:
            assert re.match(pattern, ex) is None, f"Pattern should exclude {ex}"

        # Should match other values
        for inc in should_include:
            assert re.match(pattern, inc) is not None, f"Pattern should include {inc}"

    def test_add_filter_overwrites_duplicate_keys(self):
        pf = PostFilter()
        pf.add_filter("subject", "01", None)
        assert pf.inclusions["subject"] == ["01"]
        pf.add_filter("subject", "02", None)
        assert pf.inclusions["subject"] == ["02"]

    def test_multiple_filters_in_postfilter(self):
        pf = PostFilter()
        pf.add_filter("subject", "01", None)
        pf.add_filter("session", None, "s1")
        assert "subject" in pf.inclusions
        assert "session" in pf.exclusions


class TestUnifiedFilterUnit:
    """Unit tests for UnifiedFilter class."""

    @pytest.mark.parametrize("attr", ["prefilters", "get"])
    def test_regex_search_removed_from_prefilters(self, attr: str):
        """Test that regex_search is removed from prefilters."""
        comp = {
            "filters": {"subject": "01", "regex_search": True},
            "wildcards": ["subject"],
        }
        uf = UnifiedFilter(component=comp, postfilters=PostFilter())
        compiled = getattr(uf, attr)

        assert "regex_search" not in compiled
        if attr == "prefilters":
            assert compiled == {"subject": "01"}
        else:
            assert "subject" in compiled

    @pytest.mark.parametrize(
        ("inclusion", "postfilter_key", "postfilter_val", "expected_in_get"),
        [
            (True, "subject", "01", True),
            (True, "session", "1", True),
            (False, "subject", "01", False),
            (False, "session", "1", False),
        ],
    )
    def test_postfilters_incorporated_into_unified_filter(
        self,
        inclusion: bool,
        postfilter_key: str,
        postfilter_val: str,
        expected_in_get: bool,
    ):
        pf = PostFilter()
        if inclusion:
            pf.add_filter(postfilter_key, postfilter_val, None)
        else:
            pf.add_filter(postfilter_key, None, postfilter_val)

        uf = UnifiedFilter(
            component={
                "filters": {"other": "value"},
                "wildcards": [postfilter_key],
            },
            postfilters=pf,
        )

        if expected_in_get:
            assert postfilter_key in uf.get
        else:
            assert postfilter_key in uf.post_exclusions

    def test_empty_postfilter_inclusion_turned_into_pybids_query(self):
        """Test that empty inclusion postfilter becomes Query.ANY."""
        pf = PostFilter()
        pf.add_filter("subject", [], None)

        uf = UnifiedFilter(
            component={"wildcards": ["subject"]},
            postfilters=pf,
        )

        assert uf.get["subject"][0] is Query.ANY

    @pytest.mark.parametrize(
        "inclusion",
        [True, False],
    )
    def test_postfilters_skipped_when_not_wildcards(self, inclusion: bool):
        """Test that postfilters are skipped when entity not in wildcards."""
        pf = PostFilter()
        if inclusion:
            pf.add_filter("subject", "01", None)
            pf.add_filter("session", "1", None)
        else:
            pf.add_filter("subject", None, "01")
            pf.add_filter("session", None, "1")

        # Only subject in wildcards
        uf = UnifiedFilter(
            component={"wildcards": ["subject"]},
            postfilters=pf,
        )

        if inclusion:
            assert "subject" in uf.get
            assert "session" not in uf.get
        else:
            assert "subject" in uf.post_exclusions
            assert "session" not in uf.post_exclusions

    @pytest.mark.parametrize(
        "inclusion",
        [True, False],
    )
    def test_post_exclusions_skips_entities_that_are_filters(self, inclusion: bool):
        """Test that postfilters skip entities already in prefilters."""
        pf = PostFilter()
        if inclusion:
            pf.add_filter("subject", "01", None)
            pf.add_filter("session", "1", None)
        else:
            pf.add_filter("subject", None, "01")
            pf.add_filter("session", None, "1")

        uf = UnifiedFilter(
            component={
                "wildcards": ["subject", "session"],
                "filters": {"subject": "02"},
            },
            postfilters=pf,
        )

        # Subject is in prefilters, so should use prefilter value
        if inclusion:
            assert uf.get["subject"] == ["02"]
        else:
            # Exclusions don't override prefilters
            assert "subject" not in uf.post_exclusions

    @pytest.mark.parametrize(
        ("attr", "filter_dict", "expected_keys"),
        [
            ("get", {"subject": "01", "session": {"search": "ses.*"}}, ["subject"]),
            ("search", {"subject": "01", "session": {"search": "ses.*"}}, ["session"]),
            ("get", {"sub": "01", "ses": {"match": "ses-01"}}, ["sub"]),
            ("search", {"sub": "01", "ses": {"match": "ses-01"}}, ["ses"]),
        ],
    )
    def test_get_and_search_return_appropriate_methods(
        self,
        attr: str,
        filter_dict: dict[str, str | dict[str, str]],
        expected_keys: list[str],
    ):
        """Test that get and search return correct filter methods."""
        uf = UnifiedFilter.from_filter_dict(filter_dict)  # type: ignore[arg-type]
        compiled = getattr(uf, attr)

        for key in expected_keys:
            assert key in compiled

    @pytest.mark.parametrize(
        ("empty_list", "expected"),
        [
            (True, True),
            (False, False),
        ],
    )
    def test_has_empty_prefilter_detects_empty_list(
        self,
        empty_list: bool,
        expected: bool,
    ):
        """Test detection of empty prefilter lists."""
        filters = {"subject": "01", "session": "1"}
        if empty_list:
            filters["run"] = []
        uf = UnifiedFilter.from_filter_dict(filters)
        assert uf.has_empty_prefilter == expected

    @pytest.mark.parametrize(
        ("empty_list", "expected"),
        [
            (True, True),
            (False, False),
        ],
    )
    def test_has_empty_postfilter_detects_empty_list(
        self,
        empty_list: bool,
        expected: bool,
    ):
        """Test detection of empty postfilter lists."""
        pf = PostFilter()
        pf.add_filter("subject", "01", None)
        if empty_list:
            pf.add_filter("session", [], None)
        uf = UnifiedFilter(
            component={"wildcards": ["subject", "session"]},
            postfilters=pf,
        )
        assert uf.has_empty_postfilter == expected

    @pytest.mark.parametrize(
        "value",
        [
            True,
            False,
            [True],
            [False],
            [True, "text"],
        ],
    )
    def test_without_bools_raises_when_bools_present(self, value: bool | list):
        """Test that without_bools raises when boolean filters present."""
        uf = UnifiedFilter.from_filter_dict({"subject": value})
        with pytest.raises(
            ValueError,
            match="Boolean filters in items with custom paths are not supported",
        ):
            _ = uf.without_bools

    @pytest.mark.parametrize(
        "value",
        [
            "01",
            ["01", "02"],
        ],
    )
    def test_without_bools_returns_when_no_bools(self, value: str | list[str]):
        """Test that without_bools returns mapping when no bools."""
        uf = UnifiedFilter.from_filter_dict({"subject": value})
        result = uf.without_bools
        assert "subject" in result


class TestCompileFiltersUnit:
    """Unit tests for _compile_filters function."""

    @pytest.mark.parametrize(
        ("method_a", "method_b"),
        [
            ("get", "match"),
            ("get", "search"),
            ("match", "search"),
        ],
    )
    def test_multiple_filtering_methods_raises_error(
        self,
        method_a: str,
        method_b: str,
    ):
        """Test that multiple filtering methods raise error."""
        filt = {method_a: "_", method_b: "_"}
        with pytest.raises(_TooManyKeysError):
            _compile_filters({"subject": filt}, with_regex=False)  # type: ignore[arg-type]

    def test_dict_with_no_filtering_methods_raises_error(self):
        """Test that empty dict as filter raises error."""
        with pytest.raises(_TooFewKeysError):
            _compile_filters({"subject": {}}, with_regex=False)  # type: ignore[arg-type]

    @pytest.mark.parametrize(
        "method",
        [
            "invalid",
            "regex",
            "filter",
            "query",
        ],
    )
    def test_invalid_filtering_method_raises_error(self, method: str):
        """Test that invalid filtering method raises error."""
        with pytest.raises(_InvalidKeyError):
            _compile_filters({"subject": {method: ""}}, with_regex=False)  # type: ignore[arg-type]

    @pytest.mark.parametrize(
        ("entity", "val"),
        [
            ("subject", "01"),
            ("session", "1"),
            ("run", "01"),
        ],
    )
    def test_string_filter_wrapped_with_list(self, entity: str, val: str):
        """Test that string filters are wrapped in a list."""
        res = _compile_filters({entity: val}, with_regex=False)
        assert isinstance(res, dict)
        assert entity in res
        assert isinstance(res[entity], list)
        assert res[entity][0] == val

    @pytest.mark.parametrize(
        ("kind", "with_regex", "val", "expected"),
        [
            ("simple", False, True, Query.ANY),
            ("simple", False, False, Query.NONE),
            ("list", False, True, Query.ANY),
            ("list", False, False, Query.NONE),
            ("get", False, True, Query.ANY),
            ("get", False, False, Query.NONE),
            ("match", True, True, Query.ANY),
            ("match", True, False, Query.NONE),
            ("search", True, True, Query.ANY),
            ("search", True, False, Query.NONE),
        ],
    )
    def test_boolean_filters_converted_to_pybids_queries(
        self,
        kind: str,
        with_regex: bool,
        val: bool,
        expected: object,
    ):
        """Test that boolean filters convert to Query.ANY/NONE."""
        entity = "subject"
        if kind == "simple":
            filt = {entity: val}
        elif kind == "list":
            filt = {entity: [val]}
        else:
            filt = {entity: {kind: val}}
        res = _compile_filters(filt, with_regex=with_regex)  # type: ignore[arg-type]
        assert res[entity][0] is expected

    @pytest.mark.parametrize(
        "text",
        [
            "test",
            "sub-01",
            "ses-1",
            "run-01",
        ],
    )
    def test_match_filters_converted_to_valid_regex(self, text: str):
        """Test that match filters produce anchored regex."""
        res = _compile_filters({"sub": {"match": re.escape(text)}}, with_regex=True)
        pattern = res["sub"][0]
        assert isinstance(pattern, str)
        assert re.match(pattern, text) is not None
        assert re.match(pattern, text + "_") is None

    @pytest.mark.parametrize(
        ("method", "text"),
        [
            ("match", "test"),
            ("match", "subject"),
            ("search", "test"),
            ("search", "session"),
        ],
    )
    def test_regex_filters_case_sensitive(self, method: str, text: str):
        """Test that regex filters are case sensitive."""
        res = _compile_filters({"sub": {method: text}}, with_regex=True)  # type: ignore[arg-type]
        pattern = res["sub"][0]
        assert isinstance(pattern, str)
        assert getattr(re, method)(pattern, text.upper()) is None

    @pytest.mark.parametrize(
        "text",
        [
            "test",
            "sub-01",
            "ses-1",
        ],
    )
    def test_search_filters_converted_to_valid_regex(self, text: str):
        """Test that search filters produce regex for substring matching."""
        res = _compile_filters({"sub": {"search": re.escape(text)}}, with_regex=True)
        pattern = res["sub"][0]
        assert isinstance(pattern, str)
        assert re.search(pattern, "_" + text + "_") is not None

    @pytest.mark.parametrize(
        "text",
        [
            "abc",
            "test123",
        ],
    )
    def test_search_filters_only_match_query(self, text: str):
        """Test that search filters only match the query string."""
        reversed_text = text[::-1]
        if text == reversed_text:
            pytest.skip("Text is a palindrome")
        res = _compile_filters({"sub": {"search": re.escape(text)}}, with_regex=True)
        pattern = res["sub"][0]
        assert isinstance(pattern, str)
        assert re.search(pattern, reversed_text) is None

    @pytest.mark.parametrize(
        "text",
        [
            "test",
            "sub-01",
            "",
        ],
    )
    def test_get_filters_stay_as_str(self, text: str):
        """Test that explicit get method keeps string as-is."""
        res = _compile_filters({"sub": {"get": text}}, with_regex=False)
        assert res["sub"][0] == text

    @pytest.mark.parametrize(
        ("method", "with_regex"),
        [
            ("get", True),
            ("match", False),
            ("search", False),
            (None, True),
        ],
    )
    def test_filter_skipped_when_with_regex_not_matching(
        self,
        method: str | None,
        with_regex: bool,
    ):
        """Test that filters are skipped when with_regex doesn't match."""
        entity = "subject"
        filt = {entity: ""} if method is None else {entity: {method: ""}}
        res = _compile_filters(filt, with_regex=with_regex)  # type: ignore[arg-type]
        assert entity not in res
