from __future__ import annotations

import itertools as it
import re
import string
from collections.abc import Sequence
from typing import Any, cast

import more_itertools as itx
import pytest
from bids.layout import Query
from hypothesis import assume, given
from hypothesis import strategies as st

from snakebids.core._querying import (
    _VALID_FILTER_METHODS,
    PostFilter,
    UnifiedFilter,
    _compile_filters,
    _InvalidKeyError,
    _TooFewKeysError,
    _TooManyKeysError,
    get_matching_files,
)
from snakebids.exceptions import PybidsError
from snakebids.types import FilterMap, InputConfig
from tests import strategies as sb_st


def _filters(*, min_size: int = 0):
    return st.dictionaries(
        sb_st.bids_entity().map(str), sb_st.bids_value(), min_size=min_size
    )


_mixed_filters: st.SearchStrategy[dict[str, str | list[str] | dict[str, None]]] = (
    st.dictionaries(
        sb_st.bids_entity().map(str),
        (
            st.text()
            | st.lists(st.text(), min_size=1)
            | st.sampled_from(tuple(_VALID_FILTER_METHODS)).map(lambda m: {m: None})
        ),
        max_size=6,
    )
)


class TestPostFilter:
    """Integration tests for PostFilter using property-based testing.

    Note: Most unit tests are in test_querying_unit.py. These tests focus on
    property-based validation using hypothesis for edge cases.
    """

    @given(
        inclusions=st.dictionaries(st.text(), st.text() | st.lists(st.text())),
        exclusions=st.dictionaries(st.text(), st.text() | st.lists(st.text())),
    )
    def test_postfilter_equal_to_self(
        self,
        inclusions: dict[str, str | list[str]],
        exclusions: dict[str, str | list[str]],
    ):
        """Property test: Two PostFilters with same filters should be equal."""
        pf1 = PostFilter()
        pf2 = PostFilter()
        for pf in (pf1, pf2):
            for key, val in inclusions.items():
                pf.add_filter(key, inclusions=val, exclusions=None)
            for key, val in exclusions.items():
                pf.add_filter(key, inclusions=None, exclusions=val)
        assert pf1 == pf2

    @given(key=st.text(), val=st.text() | st.iterables(st.text()))
    def test_add_filter_turns_inclusions_into_list(
        self, key: str, val: list[str] | str
    ):
        """Property test: Inclusions should always be stored as lists."""
        pf = PostFilter()
        pf.add_filter(key, val, None)
        assert isinstance(pf.inclusions[key], list)

    @given(key=st.text(), vals=st.text() | st.lists(st.text(), min_size=1))
    def test_add_filter_turns_exclusions_into_regex(self, key: str, vals: list[str]):
        """Property test: Exclusions should be stored as regex strings."""
        pf = PostFilter()
        pf.add_filter(key, None, vals)
        # exclusions are stored as a list containing a regex string
        assert key in pf.exclusions
        excl = pf.exclusions[key]
        assert isinstance(excl, list)
        assert len(excl) == 1
        assert isinstance(excl[0], str)

    @given(
        exclusion=st.text() | st.lists(st.text(), min_size=1),
        test=st.text(),
    )
    def test_exclusion_regex_excludes_correct_values(
        self, exclusion: str | list[str], test: str
    ):
        """Property test: Exclusion regex should match non-excluded values."""
        excl_list = list(itx.always_iterable(exclusion))
        assume(test not in excl_list)

        pf = PostFilter()
        pf.add_filter("", None, exclusion)
        pattern = pf.exclusions[""][0]
        # normalize exclusion into a list for testing
        # should not match excluded values
        for ex in excl_list:
            assert re.match(pattern, ex) is None
        # should match other values
        assert re.match(pattern, "".join(excl_list) + "_") is not None
        assert re.match(pattern, "_" + "".join(excl_list)) is not None
        assert re.match(pattern, test) is not None


class TestUnifiedFilter:
    """Integration tests for UnifiedFilter using property-based testing.

    Note: Most unit tests are in test_querying_unit.py. These tests focus on
    property-based validation using hypothesis for edge cases.
    """

    @pytest.mark.parametrize("inclusion", [True, False])
    @given(filters=_filters(min_size=2))
    def test_postfilters_skipped_when_not_wildcards(
        self, filters: dict[str, str], inclusion: bool
    ):
        """Property test: Postfilters should be skipped for non-wildcard entities."""
        pf = PostFilter()
        for k, v in filters.items():
            if inclusion:
                pf.add_filter(k, v, None)
            else:
                pf.add_filter(k, None, v)
        filters.popitem()
        # All but withdrawn entity in wildcards
        uf = UnifiedFilter(component={"wildcards": list(filters)}, postfilters=pf)
        assert set(filters) == set(uf.get if inclusion else uf.post_exclusions)

    @pytest.mark.parametrize("attr", ["get", "search"])
    @given(filters=_mixed_filters)
    def test_get_and_search_return_appropriate_methods(
        self,
        filters: dict[str, str | list[str] | dict[str, None]],
        attr: str,
    ):
        """Property test: get and search should return filters for appropriate methods."""
        # Build UnifiedFilter from the mixed filters and check that `get` only
        # contains 'get' (including simple string/list forms) and `search` only
        # contains 'match'/'search' filters.
        uf = UnifiedFilter.from_filter_dict(filters)  # type: ignore[arg-type]
        compiled = getattr(uf, attr)

        for key, raw in filters.items():
            method = next(iter(raw.keys())) if isinstance(raw, dict) else "get"

            if attr == "get":
                expected = method == "get"
            else:
                expected = method in ("match", "search")

            assert (key in compiled) == expected

class TestCompileFilters:
    """Integration tests for _compile_filters using property-based testing.

    Note: Most unit tests are in test_querying_unit.py. These tests focus on
    property-based validation using hypothesis for edge cases.
    """

    @given(entity=sb_st.bids_entity().map(str), val=st.text())
    def test_string_filter_wrapped_with_list(self, entity: str, val: str):
        """Property test: String filters should be wrapped into a list."""
        res = _compile_filters({entity: val}, with_regex=False)
        assert isinstance(res, dict)
        assert entity in res
        assert isinstance(res[entity], list)
        assert res[entity][0] == val

    @given(txt=st.text())
    def test_match_filters_converted_to_valid_regex(self, txt: str):
        """Property test: Match filters should produce anchored regex."""
        res = _compile_filters({"": {"match": re.escape(txt)}}, with_regex=True)
        pattern = cast(str, res[""][0])
        assert re.match(pattern, txt) is not None
        assert re.match(pattern, txt + "_") is None

    @pytest.mark.parametrize("method", ["match", "search"])
    @given(txt=st.text(alphabet=string.ascii_lowercase, min_size=1))
    def test_regex_filters_case_sensitive(self, method: str, txt: str):
        """Property test: Regex filters should be case sensitive."""
        res = _compile_filters({"": {method: txt}}, with_regex=True)  # pyright: ignore[reportArgumentType]
        pattern = cast(str, res[""][0])
        assert getattr(re, method)(pattern, txt.upper()) is None

    @given(txt=st.text())
    def test_search_filters_converted_to_valid_regex(self, txt: str):
        """Property test: Search filters should match when substring present."""
        res = _compile_filters({"": {"search": re.escape(txt)}}, with_regex=True)
        pattern = cast(str, res[""][0])
        assert re.search(pattern, "_" + txt + "_") is not None

    @given(txt=st.text(min_size=2).filter(lambda t: t != t[::-1]))
    def test_search_filters_only_match_query(self, txt: str):
        """Property test: Search filters should only match the query string."""
        res = _compile_filters({"": {"search": re.escape(txt)}}, with_regex=True)
        pattern = cast(str, res[""][0])
        assert re.search(pattern, txt[::-1]) is None


class _FakeBIDSLayout:
    def __init__(self, get: list[str] | None = None, search: list[str] | None = None):
        self.get_called = False
        self.search_called = False
        self.get_return = get or []
        self.search_return = search or []

    def get(self, regex_search: bool, **filters: Sequence[str]):
        query = re.search(
            filters.popitem()[1][0] if filters else "get", "searchgeterror"
        )
        assert query
        query = query[0]
        if query == "error":
            raise AttributeError
        assert regex_search == (query == "search")
        if query == "search":
            self.search_called = True
            return self.search_return
        assert query == "get"
        self.get_called = True
        return self.get_return


class TestGetMatchingFiles:
    def _create_unified_filter(self, *, get: bool, search: bool, error: bool):
        """Create a small UnifiedFilter for tests.

        Arguments are keyword-only to make call sites explicit about intent.
        """
        filter_dict: FilterMap = {}
        if get:
            filter_dict["get"] = "error" if error else "get"
        if search:
            filter_dict["search"] = {"search": "error" if error else "search"}
        return UnifiedFilter.from_filter_dict(filter_dict)

    @pytest.mark.parametrize("get", [False, True])
    @pytest.mark.parametrize("search", [False, True])
    def test_get_and_search_called_correctly(self, get: bool, search: bool):
        """Test appropriate call of BIDSLayout.get.

        When called with filters from UnifiedFilter.get, regex_search should be False. A
        call with regex_search == False should be made every time.

        When called with filters from UnifiedFilter.search, regex_search should be True.
        A call with regex_search == True should only be made when UnifiedFilter.search
        returns a Filter.
        """
        layout = _FakeBIDSLayout()
        filters = self._create_unified_filter(get=get, search=search, error=False)
        get_matching_files(layout, filters)  # pyright: ignore[reportArgumentType]
        assert layout.get_called
        assert layout.search_called == search

    @pytest.mark.parametrize(("get", "search"), [(True, False), (False, True)])
    def test_get_and_search_catches_attribute_error(self, get: bool, search: bool):
        layout = _FakeBIDSLayout()
        filters = self._create_unified_filter(get=get, search=search, error=True)
        with pytest.raises(
            PybidsError,
            match="Pybids has encountered a problem that Snakebids cannot handle. ",
        ):
            get_matching_files(layout, filters)  # pyright: ignore[reportArgumentType]

    @given(get=st.sets(st.text()), search=st.sets(st.text()))
    def test_returns_intersection_of_get_and_search(
        self, get: set[str], search: set[str]
    ):
        files = get_matching_files(
            _FakeBIDSLayout(get, search),  # pyright: ignore[reportArgumentType]
            self._create_unified_filter(get=True, search=True, error=False),
        )
        assert set(files) == get & search

    def test_preserves_order_of_files(self):
        get = list(map(str, range(1000)))
        search = list(map(str, range(1000)[::2]))
        layout = _FakeBIDSLayout(get, search)
        filters = self._create_unified_filter(get=True, search=True, error=False)
        files = get_matching_files(layout, filters)  # pyright: ignore[reportArgumentType]
        assert files == search
