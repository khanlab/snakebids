from __future__ import annotations

import operator as op
import re
from typing import Any

import more_itertools as itx
import pytest
from hypothesis import assume, given, settings
from hypothesis import strategies as st

import snakebids.tests.strategies as sb_st
from snakebids.utils.utils import MultiSelectDict, matches_any


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

    # This test has been flaky on github actions. No reason it should be long, so
    # just disable deadline
    @settings(deadline=None)
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
