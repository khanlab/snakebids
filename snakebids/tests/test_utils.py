from __future__ import annotations

import operator as op
import re
import string
import sys
from typing import Any

import more_itertools as itx
import pytest
from hypothesis import assume, given
from hypothesis import strategies as st

import snakebids.tests.strategies as sb_st
from snakebids.utils.utils import ImmutableList, MultiSelectDict, matches_any


class TestMatchesAny:
    @given(st.text(), st.lists(st.text()))
    def test_with_eq_operator(self, item: str, match_list: list[str]):
        if item not in match_list:
            assert matches_any(item, match_list, op.eq) is False
            match_list.append(item)
        assert matches_any(item, match_list, op.eq)

    @st.composite
    def NonMatchingMatchList(draw: st.DrawFn) -> tuple[str, list[str]]:
        item = draw(st.text())
        match_list = draw(
            st.lists(
                st.text(
                    # blank match strings (e.g. ['']) will match anything, and are not
                    # considered here
                    min_size=1
                )
                .map(re.escape)
                .filter(lambda s: not re.match(s, item))
            )
        )
        return item, match_list

    @given(NonMatchingMatchList())
    def test_with_re_match(self, args: tuple[str, list[str]]):
        item, match_list = args

        if item not in match_list:
            assert matches_any(item, match_list, re.match) is False
            match_list.append(re.escape(item))
        assert matches_any(item, match_list, re.match)

    @given(st.emails())
    def test_with_email(self, email: str):
        match_list = ["non[^@]*?email[pattern]"]
        assert matches_any(email, match_list, re.match) is False
        # Email regex copied from https://www.emailregex.com/
        match_list.append(
            r'(?:[a-z0-9!#$%&\'*+/=?^_`{|}~-]+(?:\.[a-z0-9!#$%&\'*+/=?^_`{|}~-]+)*|"(?:'
            r"[\x01-\x08\x0b\x0c\x0e-\x1f\x21\x23-\x5b\x5d-\x7f]|\\[\x01-\x09\x0b\x0c"
            r'\x0e-\x7f])*")@(?:(?:[a-z0-9](?:[a-z0-9-]*[a-z0-9])?\.)+[a-z0-9](?:[a-z0-'
            r"9-]*[a-z0-9])?|\[(?:(?:25[0-5]|2[0-4][0-9]|[01]?[0-9][0-9]?)\.){3}(?:25[0"
            r"-5]|2[0-4][0-9]|[01]?[0-9][0-9]?|[a-z0-9-]*[a-z0-9]:(?:[\x01-\x08\x0b\x0c"
            r"\x0e-\x1f\x21-\x5a\x53-\x7f]|\\[\x01-\x09\x0b\x0c\x0e-\x7f])+)\])"
        )
        assert matches_any(email, match_list, re.match, re.IGNORECASE)


class TestMultiselectDict:
    def get_selectors(
        self,
        data: st.DataObject,
        dicts: MultiSelectDict[str, Any],
        use_nonexistant_keys: bool = False,
        unique: bool = False,
    ) -> tuple[str, ...]:
        if not dicts and not use_nonexistant_keys:
            return tuple()
        sampler = st.sampled_from(list(dicts))
        val_strat = (
            st.text().filter(lambda s: s not in dicts)
            if use_nonexistant_keys
            else sampler
        )
        return tuple(data.draw(st.lists(val_strat, unique=unique)))

    @given(dicts=sb_st.multiselect_dicts(st.text(), sb_st.everything()), data=st.data())
    def test_multiple_selection_returns_same_type(
        self, dicts: MultiSelectDict[str, Any], data: st.DataObject
    ):
        selectors = self.get_selectors(data, dicts)
        # `noqa`, because in this case we actually want the classes to be exactly the
        # same
        assert type(dicts) == type(dicts[selectors])  # noqa: E721

    @given(dicts=sb_st.multiselect_dicts(st.text(), sb_st.everything()), data=st.data())
    def test_selected_items_in_original(
        self, dicts: MultiSelectDict[str, Any], data: st.DataObject
    ):
        selectors = self.get_selectors(data, dicts)
        selected = dicts[selectors]
        for key in selected:
            assert selected[key] is dicts[key]

    @given(dicts=sb_st.multiselect_dicts(st.text(), sb_st.everything()), data=st.data())
    def test_all_requested_items_received(
        self, dicts: MultiSelectDict[str, Any], data: st.DataObject
    ):
        selectors = self.get_selectors(data, dicts)
        selected = dicts[selectors]
        for selector in selectors:
            assert selector in selected

    @given(dicts=sb_st.multiselect_dicts(st.text(), sb_st.everything()), data=st.data())
    def test_no_extra_items_given(
        self, dicts: MultiSelectDict[str, Any], data: st.DataObject
    ):
        selectors = self.get_selectors(data, dicts)
        selected = dicts[selectors]
        for selector in selected:
            assert selector in selectors

    @given(
        dicts=sb_st.multiselect_dicts(st.text(), sb_st.everything(), min_size=1),
        data=st.data(),
    )
    def test_single_key_gives_single_item(
        self, dicts: MultiSelectDict[str, Any], data: st.DataObject
    ):
        selector = data.draw(st.sampled_from(list(dicts)))
        assert not isinstance(dicts[selector], MultiSelectDict)

    @given(
        dicts=sb_st.multiselect_dicts(st.text(), sb_st.everything(), min_size=1),
        data=st.data(),
    )
    def test_single_missing_key_raises_error(
        self, dicts: MultiSelectDict[str, Any], data: st.DataObject
    ):
        selector = data.draw(st.text().filter(lambda s: s not in dicts))
        with pytest.raises(KeyError):
            dicts[selector]

    @given(dicts=sb_st.multiselect_dicts(st.text(), sb_st.everything()), data=st.data())
    def test_multiple_missing_key_raises_error(
        self, dicts: MultiSelectDict[str, Any], data: st.DataObject
    ):
        selectors = self.get_selectors(data, dicts, use_nonexistant_keys=True)
        assume(len(selectors))
        with pytest.raises(KeyError):
            dicts[selectors]

    @given(dicts=sb_st.multiselect_dicts(st.text(), sb_st.everything()), data=st.data())
    def test_order_of_selectors_is_preserved(
        self, dicts: MultiSelectDict[str, Any], data: st.DataObject
    ):
        selectors = self.get_selectors(data, dicts)
        # Using the itx method for uniqueness to avoid calculating unique values in the
        # test in the same way as the source code
        assert tuple(dicts[selectors]) == tuple(itx.unique_everseen(selectors))


class TestImmutableListsAreEquivalentToTuples:
    @given(st.lists(sb_st.everything()))
    def test_they_are_not_tuples(self, items: list[Any]):
        iml = ImmutableList(items)
        assert not isinstance(iml, tuple)

    @given(st.lists(sb_st.partially_ordered(), min_size=1), st.data())
    def test_contains_items(self, items: list[Any], data: st.DataObject):
        iml = ImmutableList(items)
        sample = data.draw(st.sampled_from(items))
        assert sample in iml

    @given(st.lists(sb_st.hashables()))
    def test_is_hashable(self, items: list[Any]):
        iml = ImmutableList(items)
        assert hash(iml) == hash(tuple(items))

    @given(st.lists(sb_st.everything()))
    def test_is_iterable(self, items: list[Any]):
        i = 0
        iml = ImmutableList(items)
        _iter = iter(iml)
        while True:
            try:
                val = next(_iter)
                assert val is items[i]
                i += 1
            except StopIteration:
                break
        assert i == len(items)

    @given(st.lists(sb_st.everything()))
    def test_is_reversable(self, items: list[Any]):
        i = 0
        iml = ImmutableList(items)
        _iter = reversed(iml)
        while True:
            try:
                val = next(_iter)
                i += 1
                assert val is items[-i]
            except StopIteration:
                break
        assert i == len(items)

    @given(st.lists(sb_st.everything()))
    def test_has_length(self, items: list[Any]):
        iml = ImmutableList(items)
        assert len(iml) == len(items)

    @given(st.lists(sb_st.everything(), min_size=1), st.data())
    def test_supports_getting(self, items: list[Any], data: st.DataObject):
        iml = ImmutableList(items)
        index = data.draw(st.integers(min_value=0, max_value=len(items) - 1))
        assert iml[index] is items[index]

    @given(st.lists(sb_st.everything()), st.data())
    def test_supports_slicing(self, items: list[Any], data: st.DataObject):
        iml = ImmutableList(items)
        _slice = data.draw(st.slices(len(items)))
        for i, item in enumerate(iml[_slice]):
            assert item is items[_slice][i]

    @given(st.lists(st.integers()), st.lists(st.integers()))
    def test_supports_comparisons(self, items1: list[Any], items2: list[Any]):
        iml1 = ImmutableList(items1)
        iml2 = ImmutableList(items2)
        tup1 = tuple(items1)
        tup2 = tuple(items2)
        if tup1 == tup2:
            assert iml1 == iml2
            assert iml1 == tup2
            assert tup1 == iml2
        if tup1 >= tup2:
            assert iml1 >= iml2
            assert iml1 >= tup2
            assert tup1 >= iml2
        if tup1 > tup2:
            assert iml1 > iml2
            assert iml1 > tup2
            assert tup1 > iml2
        if tup1 < tup2:
            assert iml1 < iml2
            assert iml1 < tup2
            assert tup1 < iml2
        if tup1 <= tup2:
            assert iml1 <= iml2
            assert iml1 <= tup2
            assert tup1 <= iml2

    @given(st.lists(sb_st.partially_ordered()), st.lists(sb_st.partially_ordered()))
    def test_supports_equality(self, items1: list[Any], items2: list[Any]):
        iml1 = ImmutableList(items1)
        iml2 = ImmutableList(items2)
        tup1 = tuple(items1)
        tup2 = tuple(items2)
        if tup1 == tup2:
            assert iml1 == iml2
            assert iml1 == tup2
            assert tup1 == iml2

    @given(st.lists(sb_st.partially_ordered()))
    def test_equal_items_are_equal(self, items: list[Any]):
        iml1 = ImmutableList(items)
        iml2 = ImmutableList(items)
        tup1 = tuple(items)
        tup2 = tuple(items)
        assert iml1 == iml2
        assert iml1 == tup2
        assert tup1 == iml2

    @given(st.lists(sb_st.everything()))
    def test_supports_addition(self, items: list[Any]):
        iml = ImmutableList(items)
        assert iml + iml == tuple(items) + tuple(items)

    @given(
        st.lists(sb_st.everything(), max_size=10),
        st.integers(min_value=-5, max_value=5),
    )
    def test_supports_multiplication(self, items: list[Any], value: int):
        iml = ImmutableList(items)
        assert iml * value == tuple(items) * value
        assert value * iml == value * tuple(items)
        assert value * iml == iml * value

    @given(st.lists(sb_st.partially_ordered(), min_size=1), st.data())
    def test_supports_count(self, items: list[Any], data: st.DataObject):
        iml = ImmutableList(items)
        sample = data.draw(st.sampled_from(items))
        assert iml.count(sample) == items.count(sample)

    @given(st.lists(sb_st.partially_ordered(), min_size=1), st.data())
    def test_supports_index(self, items: list[Any], data: st.DataObject):
        iml = ImmutableList(items)
        _slice = data.draw(st.slices(len(items)))
        assume(len(items[_slice]))
        assume(_slice.step is None or _slice.step > 0)
        sample = data.draw(st.sampled_from(items[_slice]))
        start = 0 if _slice.start is None else _slice.start
        stop = sys.maxsize if _slice.stop is None else _slice.stop
        assert iml.index(sample, start, stop) == items.index(sample, start, stop)

    @given(st.lists(sb_st.everything()), sb_st.everything())
    def test_cannot_be_mutated(self, items: list[Any], value: Any):
        iml = ImmutableList(items)
        with pytest.raises(TypeError):
            iml[len(iml)] = value  # type: ignore

    @given(st.lists(st.text(string.printable)))
    def test_repr(self, items: list[Any]):
        iml = ImmutableList(items)
        assert repr(iml) == f"ImmutableList({items})"

    @given(st.lists(sb_st.everything()))
    def test_bool(self, items: list[Any]):
        iml = ImmutableList(items)
        assert bool(iml) == bool(items)
