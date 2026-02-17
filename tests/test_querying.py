from __future__ import annotations

import itertools as it
import re
import string
from collections.abc import Sequence
from typing import cast

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
    def test_postfilter_initialized_empty(self):
        pf = PostFilter()
        assert pf.inclusions == {}
        assert pf.exclusions == {}

    @given(key=st.text(), val=st.text() | st.iterables(st.text()))
    def test_add_filter_turns_inclusions_into_list(
        self, key: str, val: list[str] | str
    ):
        pf = PostFilter()
        pf.add_filter(key, val, None)
        assert isinstance(pf.inclusions[key], list)

    @given(key=st.text(), val=st.text() | st.lists(st.text()))
    def test_add_filter_copies_inclusions_without_modification(
        self, key: str, val: list[str] | str
    ):
        pf = PostFilter()
        pf.add_filter(key, val, None)
        if isinstance(val, str):
            assert pf.inclusions[key] == [val]
        else:
            assert pf.inclusions[key] == list(val)

    @given(key=st.text(), vals=st.text() | st.lists(st.text(), min_size=1))
    def test_add_filter_turns_exclusions_into_regex(self, key: str, vals: list[str]):
        pf = PostFilter()
        pf.add_filter(key, None, vals)
        # exclusions are stored as a list containing a regex string
        assert key in pf.exclusions
        excl = pf.exclusions[key]
        assert isinstance(excl, list)
        assert len(excl) == 1
        assert isinstance(excl[0], str)

    def test_empty_list_of_exclusions_treated_as_none(self):
        pf = PostFilter()
        pf.add_filter("", None, iter([]))
        assert "" not in pf.exclusions

    @given(
        exclusion=st.text() | st.lists(st.text(), min_size=1),
        test=st.text(),
    )
    def test_exclusion_regex_excludes_correct_values(
        self, exclusion: str | list[str], test: str
    ):
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

    @given(key=st.text())
    def test_add_filter_overwrites_duplicate_keys(self, key: str):
        pf = PostFilter()
        pf.add_filter(key, "a", None)
        assert pf.inclusions[key] == ["a"]
        pf.add_filter(key, "b", None)
        assert pf.inclusions[key] == ["b"]

    def test_multiple_filters_in_postfilter(self):
        pf = PostFilter()
        pf.add_filter("subject", "a", None)
        pf.add_filter("session", None, "s1")
        assert "subject" in pf.inclusions

        assert "session" in pf.exclusions


class TestUnifiedFilter:
    @pytest.mark.parametrize("attr", ["prefilters", "get"])
    @given(comp=sb_st.input_configs(wildcards=False, paths=False))
    def test_regex_search_removed_from_prefilters(self, comp: InputConfig, attr: str):
        # ensure regex_search present in component filters
        filters = dict(comp.get("filters", {}), regex_search=True)
        comp["filters"] = filters
        uf = UnifiedFilter(component=comp, postfilters=PostFilter())
        compiled = getattr(uf, attr)

        assert "regex_search" in filters
        assert "regex_search" not in compiled
        # other filters should be preserved (presence for .get, exact for .prefilters)
        filters.pop("regex_search")
        if attr == "prefilters":
            assert compiled == filters
        else:
            assert set(compiled) == set(filters)

    @pytest.mark.parametrize("inclusion", [True, False])
    @given(filters=_filters())
    def test_postfilters_incorporated_into_unified_filter(
        self, filters: dict[str, str], inclusion: bool
    ):
        pf = PostFilter()
        for k, v in filters.items():
            if inclusion:
                pf.add_filter(k, v, None)
            else:
                pf.add_filter(k, None, v)
        dummy = "".join(filters.keys()) + "_"
        uf = UnifiedFilter(
            component={
                "filters": {dummy: ""},
                "wildcards": list(filters.keys()),
            },
            postfilters=pf,
        )
        if inclusion:
            assert set(uf.get) == set(filters) | {dummy}
        else:
            assert set(uf.post_exclusions) == set(filters)

    def test_empty_postfilter_inclusion_turned_into_pybids_query(self):
        # An empty inclusion postfilter should be turned into Query.ANY in uf.get
        pf = PostFilter()
        pf.add_filter("", [], None)

        # component should list the entity as a wildcard so postfilter applies
        uf = UnifiedFilter(
            component={"wildcards": [""]},
            postfilters=pf,
        )

        assert uf.get[""][0] is Query.ANY

    @pytest.mark.parametrize("inclusion", [True, False])
    @given(filters=_filters(min_size=2))
    def test_postfilters_skipped_when_not_wildcards(
        self, filters: dict[str, str], inclusion: bool
    ):
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

    @pytest.mark.parametrize("inclusion", [True, False])
    @given(filters=_filters(min_size=2))
    def test_post_exclusions_skips_entities_that_are_filters(
        self, filters: dict[str, str], inclusion: bool
    ):
        pf = PostFilter()
        for k, v in filters.items():
            if inclusion:
                pf.add_filter(k, v, None)
            else:
                pf.add_filter(k, None, v)
        wildcards = list(filters)
        entity, value = filters.popitem()
        uf = UnifiedFilter(
            component={"wildcards": wildcards, "filters": {entity: value * 2}},
            postfilters=pf,
        )
        # since entity is present in prefilters, it should be skipped
        if inclusion:
            assert uf.get[entity] == [value * 2]
        else:
            assert entity not in uf.post_exclusions

    @pytest.mark.parametrize("attr", ["get", "search"])
    @given(filters=_mixed_filters)
    def test_get_and_search_return_appropriate_methods(
        self,
        filters: dict[str, str | list[str] | dict[str, None]],
        attr: str,
    ):
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

    @pytest.mark.parametrize("empty_list", [True, False])
    @given(
        filters=st.dictionaries(
            sb_st.bids_entity().map(str), st.text() | st.lists(st.text(), min_size=1)
        )
    )
    def test_has_empty_prefilter_detects_empty_list(
        self, filters: dict[str, str | list[str]], empty_list: bool
    ):
        if empty_list:
            filters["".join(filters)] = []
        uf = UnifiedFilter.from_filter_dict(filters)
        assert uf.has_empty_prefilter == empty_list

    @pytest.mark.parametrize("empty_list", [True, False])
    @given(
        filters=st.dictionaries(
            sb_st.bids_entity().map(str), st.text() | st.lists(st.text(), min_size=1)
        )
    )
    def test_has_empty_postfilter_detects_empty_list(
        self, filters: dict[str, str | list[str]], empty_list: bool
    ):
        pf = PostFilter()
        for k, v in filters.items():
            pf.add_filter(k, v, None)
        dummy = "".join(filters)
        if empty_list:
            pf.add_filter(dummy, [], None)
        uf = UnifiedFilter(component={"wildcards": [*filters, dummy]}, postfilters=pf)
        assert uf.has_empty_postfilter == empty_list

    # without_bools behavior tested below with parameterized cases
    @given(
        value=st.booleans()
        | st.lists(st.booleans(), min_size=1)
        | st.tuples(st.booleans(), st.text()).map(list)
    )
    def test_without_bools_raises_when_bools_present(self, value: bool):
        """Accessing without_bools must raise when boolean filters are present."""
        uf = UnifiedFilter.from_filter_dict({"": value})
        with pytest.raises(
            ValueError,
            match=("Boolean filters in items with custom paths are not supported"),
        ):
            _ = uf.without_bools

    @given(value=st.text() | st.lists(st.text()))
    def test_without_bools_returns_when_no_bools(self, value: str | list[str]):
        """without_bools should return the mapping when no boolean filters exist."""
        uf = UnifiedFilter.from_filter_dict({"": value})
        result = uf.without_bools
        assert "" in result


class TestCompileFilters:
    @pytest.mark.parametrize(
        ("method_a", "method_b"),
        list(it.combinations(tuple(_VALID_FILTER_METHODS), 2)),
    )
    def test_multiple_filtering_methods_raises_error(
        self, method_a: str, method_b: str
    ):
        # dict with two different filtering methods should raise a _TooManyKeysError
        filt = {method_a: "_", method_b: "_"}
        with pytest.raises(_TooManyKeysError):
            _compile_filters({"": filt}, with_regex=False)  # type: ignore[arg-type]

    def test_dict_with_no_filtering_methods_raises_error(self):
        # empty dict as a filter value should raise a _TooFewKeysError
        with pytest.raises(_TooFewKeysError):
            _compile_filters({"": {}}, with_regex=False)  # type: ignore[arg-type]

    @given(method=st.text().filter(lambda t: t not in _VALID_FILTER_METHODS))
    def test_invalid_filtering_method_raises_error(self, method: str):
        # unknown filtering key should raise a _InvalidKeyError
        with pytest.raises(_InvalidKeyError):
            _compile_filters({"": {method: ""}}, with_regex=False)  # type: ignore[arg-type]

    @given(entity=sb_st.bids_entity().map(str), val=st.text())
    def test_string_filter_wrapped_with_list(self, entity: str, val: str):
        # plain string filters should be wrapped into a single-element list
        res = _compile_filters({entity: val}, with_regex=False)
        assert isinstance(res, dict)
        assert entity in res
        assert isinstance(res[entity], list)
        assert res[entity][0] == val

    @pytest.mark.parametrize(
        ("kind", "with_regex"),
        [
            ("simple", False),
            ("list", False),
            ("get", False),
            ("match", True),
            ("search", True),
        ],
    )
    @pytest.mark.parametrize(
        ("val", "expected"),
        [
            (True, Query.ANY),
            (False, Query.NONE),
        ],
    )
    @given(entity=sb_st.bids_entity().map(str))
    def test_boolean_filters_converted_to_pybids_queries(
        self, entity: str, kind: str, with_regex: bool, val: bool, expected: object
    ):
        """Boolean filters (True/False) should convert to Query.ANY / Query.NONE

        Test three shapes for each boolean: simple boolean, boolean in a list, and
        boolean nested under each of the valid methods ('get','match','search').
        """
        if kind == "simple":
            filt = {entity: val}
        elif kind == "list":
            filt = {entity: [val]}
        else:
            filt = {entity: {kind: val}}
        res = _compile_filters(filt, with_regex=with_regex)  # type: ignore[arg-type]
        assert res[entity][0] is expected

    @given(txt=st.text())
    def test_match_filters_converted_to_valid_regex(self, txt: str):
        # match should produce anchored regex that matches exactly
        res = _compile_filters({"": {"match": re.escape(txt)}}, with_regex=True)
        pattern = cast(str, res[""][0])
        assert re.match(pattern, txt) is not None
        assert re.match(pattern, txt + "_") is None

    @pytest.mark.parametrize("method", ["match", "search"])
    @given(txt=st.text(alphabet=string.ascii_lowercase, min_size=1))
    def test_regex_filters_case_sensitive(self, method: str, txt: str):
        res = _compile_filters({"": {method: txt}}, with_regex=True)  # pyright: ignore[reportArgumentType]
        pattern = cast(str, res[""][0])
        assert getattr(re, method)(pattern, txt.upper()) is None

    @given(txt=st.text())
    def test_search_filters_converted_to_valid_regex(self, txt: str):
        # search should produce regex that matches when substring present
        res = _compile_filters({"": {"search": re.escape(txt)}}, with_regex=True)
        pattern = cast(str, res[""][0])
        assert re.search(pattern, "_" + txt + "_") is not None

    @given(txt=st.text(min_size=2).filter(lambda t: t != t[::-1]))
    def test_search_filters_only_match_query(self, txt: str):
        # search should produce regex that matches when substring present
        res = _compile_filters({"": {"search": re.escape(txt)}}, with_regex=True)
        pattern = cast(str, res[""][0])
        assert re.search(pattern, txt[::-1]) is None

    @given(txt=st.text())
    def test_get_filters_stay_as_str(self, txt: str):
        # explicit get method should leave the string as-is
        res = _compile_filters({"": {"get": txt}}, with_regex=False)
        assert cast(str, res[""][0]) == txt

    @pytest.mark.parametrize(
        ("method", "with_regex"),
        [
            ("get", True),
            ("match", False),
            ("search", False),
            (None, True),
        ],
    )
    @given(entity=sb_st.bids_entity().map(str))
    def test_filter_skipped_when_with_regex_not_matching(
        self, entity: str, method: str | None, with_regex: bool
    ):
        filt = {entity: ""} if method is None else {entity: {method: ""}}

        res = _compile_filters(filt, with_regex=with_regex)  # type: ignore[arg-type]

        assert entity not in res


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
    def test_returns_union_of_get_and_search(self, get: set[str], search: set[str]):
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
