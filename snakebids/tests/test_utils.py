import operator as op
import re
from typing import List

from hypothesis import given
from hypothesis import strategies as st

from snakebids.utils.utils import matches_any


@given(st.text(), st.lists(st.text()))
def test_matches_any_with_eq_operator(item: str, match_list: List[str]):
    if item not in match_list:
        assert matches_any(item, match_list, op.eq) is False
        match_list.append(item)
    assert matches_any(item, match_list, op.eq)


@given(st.text(), st.lists(st.text(min_size=1).map(re.escape)))
def test_matches_any_with_re_match(item: str, match_list: List[str]):
    # blank match strings (e.g. ['']) will match anything, and are not considered here
    if item not in match_list:
        assert matches_any(item, match_list, re.match) is False
        match_list.append(re.escape(item))
    assert matches_any(item, match_list, re.match)


@given(st.emails())
def test_matches_any_with_email(email: str):
    match_list = ["non[^@]*?email[pattern]"]
    assert matches_any(email, match_list, re.match) is False
    # Email regex copied from https://www.emailregex.com/
    match_list.append(
        r'(?:[a-z0-9!#$%&\'*+/=?^_`{|}~-]+(?:\.[a-z0-9!#$%&\'*+/=?^_`{|}~-]+)*|"(?:['
        r"\x01-\x08\x0b\x0c\x0e-\x1f\x21\x23-\x5b\x5d-\x7f]|\\[\x01-\x09\x0b\x0c\x0e-"
        r'\x7f])*")@(?:(?:[a-z0-9](?:[a-z0-9-]*[a-z0-9])?\.)+[a-z0-9](?:[a-z0-9-]*'
        r"[a-z0-9])?|\[(?:(?:25[0-5]|2[0-4][0-9]|[01]?[0-9][0-9]?)\.){3}(?:25[0-5]|2"
        r"[0-4][0-9]|[01]?[0-9][0-9]?|[a-z0-9-]*[a-z0-9]:(?:[\x01-\x08\x0b\x0c\x0e-"
        r"\x1f\x21-\x5a\x53-\x7f]|\\[\x01-\x09\x0b\x0c\x0e-\x7f])+)\])"
    )
    assert matches_any(email, match_list, re.match, re.IGNORECASE)
