from __future__ import annotations

import copy
from typing import Any

import attrs
import more_itertools as itx
import pytest
from hypothesis import assume, given

from snakebids.core._table import BidsTable
from snakebids.core.datasets import BidsComponent
from snakebids.tests import strategies as sb_st
from snakebids.tests.helpers import get_bids_path
from snakebids.utils.utils import BidsEntity


@given(sb_st.bids_tables(min_entities=2))
def test_zip_lists_must_be_same_length(table: BidsTable):
    zip_lists = table.to_dict()
    itx.first(zip_lists.values()).append("foo")
    with pytest.raises(
        ValueError, match="each entity must have the same number of values"
    ):
        BidsComponent(name="foo", path=get_bids_path(zip_lists), table=zip_lists)


class TestEq:
    @given(sb_st.bids_tables(), sb_st.everything_except(BidsTable))
    def test_other_types_are_unequal(self, table: BidsTable, other: Any):
        assert table != other

    @given(sb_st.bids_tables(min_entities=0, min_values=0))
    def test_copied_object_is_equal(self, table: BidsTable):
        other = copy.deepcopy(table)
        assert table == other

    @given(sb_st.bids_tables(min_entities=2, min_values=0))
    def test_wildcard_order_is_irrelevant(self, table: BidsTable):
        other = copy.deepcopy(table)
        reordered = BidsTable(
            wildcards=reversed(other.wildcards),
            entries=[tuple(reversed(entry)) for entry in other.entries],
        )
        assert table == reordered

    @given(
        table=sb_st.bids_tables(min_entities=1, min_values=0),
        wildcard=sb_st.bids_entity(),
    )
    def test_wildcards_must_be_the_same(self, table: BidsTable, wildcard: BidsEntity):
        assume(wildcard.wildcard not in table.wildcards)
        other = copy.deepcopy(table)
        reordered = attrs.evolve(other, wildcards=[wildcard, *other.wildcards[1:]])
        assert table != reordered

    @given(sb_st.bids_tables())
    def test_mutation_makes_unequal(self, table: BidsTable):
        cp = copy.deepcopy(table)
        other = attrs.evolve(
            cp,
            entries=[
                ("0" + "".join(cp.entries[0]), *cp.entries[0][1:]),
                *cp.entries[1:],
            ],
        )
        assert table != other

    @given(sb_st.bids_tables())
    def test_extra_entry_makes_unequal(self, table: BidsTable):
        cp = copy.deepcopy(table)
        other = attrs.evolve(cp, entries=[cp.entries[0], *cp.entries])
        assert table != other

    @given(sb_st.bids_tables())
    def test_missing_entry_makes_unequal(self, table: BidsTable):
        cp = copy.deepcopy(table)
        other = attrs.evolve(cp, entries=cp.entries[1:])
        assert table != other
