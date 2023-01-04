import copy
import itertools as it
from typing import Any, Dict, List

import more_itertools as itx
import pytest
from hypothesis import assume, given
from hypothesis import strategies as st

from snakebids.core.datasets import BidsComponent, BidsDataset
from snakebids.tests import strategies as sb_st
from snakebids.tests.helpers import debug, get_bids_path, get_zip_list, setify
from snakebids.utils import sb_itertools as sb_it
from snakebids.utils.utils import BidsEntity


class TestBidsComponentAliases:
    @given(sb_st.bids_components())
    def test_bids_component_aliases_are_correctly_set(self, component: BidsComponent):
        assert component.input_path is component.path
        assert component.input_zip_lists is component.zip_lists
        assert component.input_lists is component.input_lists
        assert component.input_wildcards is component.wildcards

    @given(sb_st.bids_components())
    def test_bids_dataset_aliases_are_correctly_set(self, component: BidsComponent):
        dataset = BidsDataset.from_iterable([component])
        assert dataset.input_path == dataset.path
        assert dataset.input_zip_lists == dataset.zip_lists
        assert dataset.input_lists == dataset.input_lists
        assert dataset.input_wildcards == dataset.wildcards


class TestBidsComponentValidation:
    @given(sb_st.input_zip_lists().filter(lambda v: len(v) > 1))
    def test_zip_lists_must_be_same_length(self, zip_lists: Dict[str, List[str]]):
        itx.first(zip_lists.values()).append("foo")
        with pytest.raises(ValueError) as err:
            BidsComponent("foo", get_bids_path(zip_lists), zip_lists)
        assert err.value.args[0] == "zip_lists must all be of equal length"

    @given(sb_st.input_zip_lists(), sb_st.bids_entity())
    def test_path_cannot_have_extra_entities(
        self, zip_lists: Dict[str, List[str]], entity: BidsEntity
    ):
        assume(entity.wildcard not in zip_lists)
        path = get_bids_path(it.chain(zip_lists, [entity.entity]))
        with pytest.raises(ValueError) as err:
            BidsComponent("foo", path, zip_lists)
        assert (
            "zip_lists entries must match the wildcards in input_path"
            in err.value.args[0]
        )

    @given(sb_st.input_zip_lists().filter(lambda v: len(v) > 1))
    def test_path_cannot_have_missing_entities(self, zip_lists: Dict[str, List[str]]):
        path = get_bids_path(sb_it.drop(1, zip_lists))
        with pytest.raises(ValueError) as err:
            BidsComponent("foo", path, zip_lists)
        assert (
            "zip_lists entries must match the wildcards in input_path"
            in err.value.args[0]
        )


class TestBidsComponentEq:
    @given(sb_st.bids_components(), sb_st.everything_except(BidsComponent))
    def test_other_types_are_unequal(self, input: BidsComponent, other: Any):
        assert input != other

    def test_empty_BidsInput_are_equal(self):
        assert BidsComponent("", "", {}) == BidsComponent("", "", {})
        assert BidsComponent("", "{foo}{bar}", {"foo": [], "bar": []}) == BidsComponent(
            "", "{foo}{bar}", {"foo": [], "bar": []}
        )

    @given(sb_st.bids_components())
    def test_copies_are_equal(self, input: BidsComponent):
        cp = copy.deepcopy(input)
        assert cp == input

    @given(sb_st.bids_components())
    def test_mutation_makes_unequal(self, input: BidsComponent):
        cp = copy.deepcopy(input)
        itx.first(cp.zip_lists.values())[0] += "foo"
        assert cp != input

    @given(sb_st.bids_components(), st.data())
    def test_extra_entities_makes_unequal(
        self, input: BidsComponent, data: st.DataObject
    ):
        cp = copy.deepcopy(input)
        new_entity = data.draw(
            sb_st.bids_value().filter(lambda s: s not in input.zip_lists)
        )
        cp.zip_lists[new_entity] = []
        itx.first(cp.zip_lists.values())[0] += "foo"
        assert cp != input

    @given(sb_st.bids_components())
    def test_order_doesnt_affect_equality(self, input: BidsComponent):
        cp = copy.deepcopy(input)
        for l in cp.zip_lists:
            cp.zip_lists[l].reverse()
        assert cp == input

    @given(sb_st.bids_components())
    def test_paths_must_be_identical(self, input: BidsComponent):
        cp = BidsComponent(input.input_name, input.input_path + "foo", input.zip_lists)
        assert cp != input


class TestBidsComponentProperties:
    @given(st.data(), st.integers(min_value=1, max_value=2))
    def test_input_lists_derives_from_zip_lists(
        self, data: st.DataObject, min_size: int
    ):
        input_lists = data.draw(sb_st.bids_input_lists(min_size, max_size=5))

        # Due to the product, we can delete some of the combinations and still
        # regenerate our input_lists
        combs = list(it.product(*input_lists.values()))[min_size - 1 :]
        zip_lists = get_zip_list(input_lists, combs)
        path = get_bids_path(zip_lists)

        assert setify(
            BidsComponent(name="foo", path=path, zip_lists=zip_lists).entities
        ) == setify(input_lists)

    @given(
        st.dictionaries(
            sb_st.bids_entity().map(lambda e: e.wildcard),
            sb_st.bids_value("[^.]*"),
            min_size=1,
        ).filter(lambda v: list(v) != ["datatype"])
    )
    def test_input_wildcards_derives_from_zip_lists(
        self,
        bids_entities: Dict[str, str],
    ):
        zip_lists = {entity: [val] for entity, val in bids_entities.items()}
        bids_input = BidsComponent(
            name="foo",
            path=get_bids_path(zip_lists),
            zip_lists=zip_lists,
        )

        wildstr = ".".join(bids_input.input_wildcards.values())
        first = wildstr.format(**bids_input.input_wildcards)
        second = first.format(**bids_entities)
        assert set(second.split(".")) == set(bids_entities.values())
