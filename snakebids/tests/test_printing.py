# ruff: noqa: PLR2004
from __future__ import annotations

import pyparsing as pp
from hypothesis import assume, given
from hypothesis import strategies as st

import snakebids.tests.strategies as sb_st
from snakebids.core.datasets import BidsComponent, BidsDataset
from snakebids.io.printing import format_zip_lists
from snakebids.types import ZipList


def zip_list_parser() -> pp.ParserElement:
    key = pp.quoted_string().set_results_name("key") + pp.Suppress(":")
    elided_list = (
        pp.Opt(pp.delimited_list(pp.quoted_string()))("left")
        + pp.Opt(",").suppress()
        + pp.Opt("...")("ellipse")
        + pp.Opt(pp.delimited_list(pp.quoted_string()))("right")
    )
    row = key + pp.Suppress("[") + elided_list + pp.Suppress("],")
    return pp.Suppress("{") + pp.Group(row)[1, ...] + pp.Suppress("}")


@given(zip_list=sb_st.zip_lists(max_entities=1, restrict_patterns=True))
def test_ellipses_appears_when_maxwidth_too_short(zip_list: ZipList):
    width = len(format_zip_lists(zip_list, tabstop=0).splitlines()[1])
    parsed = zip_list_parser().parse_string(
        format_zip_lists(zip_list, width - 1, tabstop=0)
    )
    assert "ellipse" in parsed[0]


@given(
    zip_list=sb_st.zip_lists(
        min_values=1, max_values=4, max_entities=4, restrict_patterns=True
    ),
    width=st.integers(min_value=10, max_value=200),
)
def test_values_balanced_around_elision_correctly(zip_list: ZipList, width: int):
    parsed: pp.ParseResults = zip_list_parser().parse_string(
        format_zip_lists(zip_list, max_width=width, tabstop=0)
    )
    assert parsed
    assert parsed[0]
    nleft = len(parsed[0].get("left", []))  # type: ignore
    nright = len(parsed[0].get("right", []))  # type: ignore
    if "ellipse" in parsed[0]:
        assert nleft == nright or nleft == nright + 1
        assert nleft + nright == len(parsed[0]) - 2  # type: ignore
    else:
        assert nright == 0
        # nleft contains all elements in the first row except the key
        assert nleft == len(parsed[0]) - 1  # type: ignore

    for line in parsed[1:]:  # type: ignore
        assert len(line.get("left", [])) == nleft  # type: ignore
        assert len(line.get("right", [])) == nright  # type: ignore


class TestCorrectNumberOfLinesCreated:
    @given(
        zip_list=sb_st.zip_lists(max_values=1, max_entities=6, restrict_patterns=True),
    )
    def test_in_zip_list(self, zip_list: ZipList):
        assert (
            len(format_zip_lists(zip_list, tabstop=0).splitlines()) == len(zip_list) + 2
        )

    @given(
        component=sb_st.bids_components(
            max_values=1, max_entities=6, restrict_patterns=True
        ),
    )
    def test_in_component(self, component: BidsComponent):
        assert (
            len(component.pformat(tabstop=0).splitlines())
            == len(component.entities) + 6
        )

    @given(dataset=sb_st.datasets())
    def test_in_dataset(self, dataset: BidsDataset):
        n_comps = len(dataset)
        n_entities = sum(len(comp.entities) for comp in dataset.values())
        assert (
            len(dataset.pformat(tabstop=0).splitlines()) == n_entities + 6 * n_comps + 2
        )


class TestIsValidPython:
    @given(zip_list=sb_st.zip_lists(restrict_patterns=True))
    def test_in_zip_list(self, zip_list: ZipList):
        assert eval(format_zip_lists(zip_list)) == zip_list

    @given(component=sb_st.bids_components(restrict_patterns=True))
    def test_in_component(self, component: BidsComponent):
        assert eval((component.pformat())) == component

    @given(dataset=sb_st.datasets())
    def test_in_dataset(self, dataset: BidsDataset):
        assert eval(dataset.pformat()) == dataset


# this could also be tested for components and datasets, however, in those objects the
# path and name are allowed to be longer than the width, so finding the zip_list lines
# would prove more challenging than it's worth
@given(
    zip_list=sb_st.zip_lists(max_entities=1, restrict_patterns=True),
    width=st.integers(10, 100),
    tab=st.integers(0, 10),
)
def test_line_never_longer_than_max_width(zip_list: ZipList, width: int, tab: int):
    assume(width > tab + 10)
    formatted = format_zip_lists(zip_list, width, tab)
    parsed = zip_list_parser().parse_string(formatted)
    assume("left" in parsed[0])
    assert all(list(len(line) <= width for line in formatted.splitlines()))


def get_indent_length(line: str):
    return len(line) - len(line.lstrip(" "))


class TestIndentLengthMultipleOfTabStop:
    @given(zip_list=sb_st.zip_lists(restrict_patterns=True), tabstop=st.integers(1, 10))
    def test_in_zip_list(self, zip_list: ZipList, tabstop: int):
        for line in format_zip_lists(zip_list, tabstop=tabstop).splitlines():
            assert get_indent_length(line) / tabstop in {0, 1}

    @given(
        component=sb_st.bids_components(max_values=1, restrict_patterns=True),
        tabstop=st.integers(1, 10),
    )
    def test_in_component(self, component: BidsComponent, tabstop: int):
        for line in component.pformat(tabstop=tabstop).splitlines():
            assert get_indent_length(line) / tabstop in {0, 1, 2}

    @given(dataset=sb_st.datasets(), tabstop=st.integers(1, 10))
    def test_in_dataset(self, dataset: BidsDataset, tabstop: int):
        for line in dataset.pformat(tabstop=tabstop).splitlines():
            assert get_indent_length(line) / tabstop in {0, 1, 2, 3}


class TestMultipleLevelsOfIndentationUsed:
    @given(
        component=sb_st.bids_components(max_values=1, restrict_patterns=True),
        tabstop=st.integers(1, 10),
    )
    def test_in_component(self, component: BidsComponent, tabstop: int):
        indents: set[int] = set()
        for line in component.pformat(tabstop=tabstop).splitlines():
            indents.add(get_indent_length(line))
        assert len(indents) == 3

    @given(dataset=sb_st.datasets(), tabstop=st.integers(1, 10))
    def test_in_dataset(self, dataset: BidsDataset, tabstop: int):
        indents: set[int] = set()
        for line in dataset.pformat(tabstop=tabstop).splitlines():
            indents.add(get_indent_length(line))
        assert len(indents) == 4
