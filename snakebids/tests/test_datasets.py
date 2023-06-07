from __future__ import annotations

import copy
import itertools as it
import re
import string
import warnings
from pathlib import Path
from typing import Any

import more_itertools as itx
import pytest
from hypothesis import assume, given
from hypothesis import strategies as st
from snakemake.exceptions import WildcardError

from snakebids.core.construct_bids import bids
from snakebids.core.datasets import (
    BidsComponent,
    BidsComponentRow,
    BidsDataset,
    BidsPartialComponent,
)
from snakebids.exceptions import DuplicateComponentError
from snakebids.tests import strategies as sb_st
from snakebids.tests.helpers import expand_zip_list, get_bids_path, get_zip_list, setify
from snakebids.types import Expandable, ZipList
from snakebids.utils import sb_itertools as sb_it
from snakebids.utils.snakemake_io import glob_wildcards
from snakebids.utils.utils import BidsEntity, get_wildcard_dict, zip_list_eq


def test_multiple_components_cannot_have_same_name():
    comp1 = BidsComponent(name="foo", path=".", zip_lists={})
    comp2 = BidsComponent(name="foo", path=".", zip_lists={})
    with pytest.raises(DuplicateComponentError):
        BidsDataset.from_iterable([comp1, comp2])


class TestBidsComponentAliases:
    @given(sb_st.bids_components())
    def test_bids_component_aliases_are_correctly_set(self, component: BidsComponent):
        assert component.input_path is component.path
        assert component.input_zip_lists is component.zip_lists
        assert component.input_lists is component.input_lists
        assert component.input_wildcards is component.wildcards

    @given(sb_st.bids_components())
    def test_bids_dataset_aliases_are_correctly_set(self, component: BidsComponent):
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", category=DeprecationWarning)
            dataset = BidsDataset.from_iterable([component])
            assert dataset.input_path == dataset.path
            assert dataset.input_zip_lists == dataset.zip_lists
            assert dataset.input_lists == dataset.entities
            assert dataset.input_wildcards == dataset.wildcards


class TestBidsComponentValidation:
    @given(sb_st.zip_lists().filter(lambda v: len(v) > 1))
    def test_zip_lists_must_be_same_length(self, zip_lists: ZipList):
        itx.first(zip_lists.values()).append("foo")
        with pytest.raises(ValueError) as err:
            BidsComponent(
                name="foo", path=get_bids_path(zip_lists), zip_lists=zip_lists
            )
        assert err.value.args[0] == "zip_lists must all be of equal length"

    @given(sb_st.zip_lists(), sb_st.bids_entity())
    def test_path_cannot_have_extra_entities(
        self, zip_lists: ZipList, entity: BidsEntity
    ):
        assume(entity.wildcard not in zip_lists)
        path = get_bids_path(it.chain(zip_lists, [entity.entity]))
        with pytest.raises(ValueError) as err:
            BidsComponent(name="foo", path=path, zip_lists=zip_lists)
        assert (
            "zip_lists entries must match the wildcards in input_path"
            in err.value.args[0]
        )

    @given(sb_st.zip_lists().filter(lambda v: len(v) > 1))
    def test_path_cannot_have_missing_entities(self, zip_lists: ZipList):
        # Snakebids strategies won't return a zip_list with just datatype, but now that
        # we've dropped an entity we need to check again
        path_entities = sb_it.drop(1, zip_lists)
        assume(set(path_entities) - {"datatype"})

        path = get_bids_path(path_entities)
        with pytest.raises(ValueError) as err:
            BidsComponent(name="foo", path=path, zip_lists=zip_lists)
        assert (
            "zip_lists entries must match the wildcards in input_path"
            in err.value.args[0]
        )


class TestBidsComponentEq:
    @given(sb_st.bids_components(), sb_st.everything_except(BidsComponent))
    def test_other_types_are_unequal(self, input: BidsComponent, other: Any):
        assert input != other

    def test_empty_BidsInput_are_equal(self):
        assert BidsComponent(name="", path="", zip_lists={}) == BidsComponent(
            name="", path="", zip_lists={}
        )
        assert BidsComponent(
            name="",
            path="{foo}{bar}",
            zip_lists={"foo": [], "bar": []},
        ) == BidsComponent(
            name="",
            path="{foo}{bar}",
            zip_lists={"foo": [], "bar": []},
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
        for list_ in cp.zip_lists:
            cp.zip_lists[list_].reverse()
        assert cp == input

    @given(sb_st.bids_components())
    def test_paths_must_be_identical(self, input: BidsComponent):
        cp = BidsComponent(
            name=input.input_name,
            path=input.input_path + "foo",
            zip_lists=input.zip_lists,
        )
        assert cp != input


class TestBidsComponentProperties:
    @given(st.data(), st.integers(min_value=1, max_value=2))
    def test_input_lists_derives_from_zip_lists(
        self, data: st.DataObject, min_size: int
    ):
        input_lists = data.draw(sb_st.bids_input_lists(min_size=min_size, max_size=5))

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
        bids_entities: dict[str, str],
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


legacy_names = [
    "input_path",
    "input_zip_lists",
    "input_lists",
    "input_wildcards",
    "subjects",
    "sessions",
    "subj_wildcards",
]


class TestBidsDatasetLegacyAccess:
    match_str = r".*generate_inputs\(\) no longer returns a dict by default.*"

    @given(dataset=sb_st.datasets(), name=st.sampled_from(legacy_names))
    def test_accessing_legacy_names_leads_to_informative_error(
        self, dataset: BidsDataset, name: str
    ):
        assume(name not in dataset)
        with pytest.raises(KeyError, match=self.match_str):
            dataset[name]

    @given(
        dataset=sb_st.datasets(),
        name=st.text(sb_st.alphanum).filter(lambda s: s not in legacy_names),
    )
    def test_accessing_ordinary_missing_names_leads_to_regular_error(
        self, dataset: BidsDataset, name: str
    ):
        assume(name not in dataset)
        with pytest.raises(KeyError) as err:
            dataset[name]

        assert not re.search(self.match_str, str(err.value))

    @given(
        dataset=sb_st.datasets_one_comp(
            names=st.shared(st.sampled_from(legacy_names), key="names")
        ),
        name=st.shared(st.sampled_from(legacy_names), key="names"),
    )
    def test_legacy_names_can_be_accessed_if_component_name(
        self, dataset: BidsDataset, name: str
    ):
        assert isinstance(dataset[name], BidsComponent)


def _get_novel_path(prefix: str, component: Expandable):
    # use the "comp-" prefix to give a constant part to the novel template,
    # otherwise the trivial template "{foo}" globs everything
    return Path(
        *map(
            lambda s: f"{prefix}-{s}",
            get_wildcard_dict(component.zip_lists).values(),
        )
    )


class TestExpandables:
    @given(
        component=sb_st.expandables(
            restrict_patterns=True,
            blacklist_entities=["extension"],
        ),
        wildcards=st.lists(
            st.text(string.ascii_letters, min_size=1, max_size=10).filter(
                lambda s: s not in sb_st.valid_entities
            ),
            min_size=1,
            max_size=5,
        ),
        data=st.data(),
    )
    def test_expand_with_extra_args_returns_all_paths(
        self, component: Expandable, wildcards: list[str], data: st.DataObject
    ):
        num_wildcards = len(wildcards)
        values = data.draw(
            st.lists(
                st.lists(
                    st.text(sb_st.alphanum, min_size=1, max_size=10),
                    min_size=1,
                    max_size=3,
                ),
                min_size=num_wildcards,
                max_size=num_wildcards,
            )
        )
        path_tpl = bids(
            **get_wildcard_dict(component.zip_lists),
            **get_wildcard_dict(wildcards),
        )
        wcard_dict = dict(zip(wildcards, values))
        zlist = expand_zip_list(component.zip_lists, wcard_dict)
        paths = component.expand(path_tpl, **wcard_dict)
        assert zip_list_eq(glob_wildcards(path_tpl, paths), zlist)

    @given(
        component=sb_st.expandables(
            restrict_patterns=True, blacklist_entities=["extension"]
        )
    )
    def test_expand_over_multiple_paths(self, component: Expandable):
        path1 = _get_novel_path("first", component)
        path2 = _get_novel_path("second", component)
        paths = component.expand([path1, path2])
        assert zip_list_eq(glob_wildcards(path1, paths), glob_wildcards(path2, paths))

    @given(
        component=sb_st.expandables(restrict_patterns=True),
        wildcard=st.text(string.ascii_letters, min_size=1, max_size=10).filter(
            lambda s: s not in sb_st.valid_entities
        ),
    )
    def test_partial_expansion(self, component: Expandable, wildcard: str):
        path_tpl = bids(
            **get_wildcard_dict(component.zip_lists), **get_wildcard_dict(wildcard)
        )
        paths = component.expand(path_tpl, allow_missing=True)
        for path in paths:
            assert re.search(r"\{.+\}", path)

    @given(
        component=sb_st.expandables(restrict_patterns=True),
        wildcard=st.text(string.ascii_letters, min_size=1, max_size=10).filter(
            lambda s: s not in sb_st.valid_entities
        ),
    )
    def test_prevent_partial_expansion(self, component: Expandable, wildcard: str):
        path_tpl = bids(
            **get_wildcard_dict(component.zip_lists), **get_wildcard_dict(wildcard)
        )
        with pytest.raises(WildcardError):
            component.expand(path_tpl)


class TestBidsComponentExpand:
    """
    `extension` is generally excluded from the generated components because its wildcard
    is immediately adjacent to the previous wildcard (e.g.
    sub-{subject}_{suffix}{extenions}), making it impossible to correctly parse using
    `glob_wildcards`. Its absence does not affect the spirit of the tests.

    We also exclude extra_entities from the generated Components for the same reason:
    ``glob_wildcards`` would glob the extra entities found in the generated path, making
    comparison impossible
    """

    @given(
        component=sb_st.bids_components(
            blacklist_entities=["extension"],
            restrict_patterns=True,
            extra_entities=False,
        )
    )
    def test_expand_with_no_args_returns_initial_paths(self, component: BidsComponent):
        paths = component.expand()
        assert zip_list_eq(glob_wildcards(component.path, paths), component.zip_lists)

    @given(
        component=sb_st.bids_components(restrict_patterns=True, extra_entities=False)
    )
    def test_not_expand_over_internal_path_when_novel_given(
        self, component: BidsComponent
    ):
        assume(set(component.wildcards) - {"datatype", "suffix", "extension"})
        novel_path = _get_novel_path("comp", component)
        paths = component.expand(novel_path)
        assert not glob_wildcards(component.path, paths)


class TestFiltering:
    def get_filter_dict(
        self,
        data: st.DataObject,
        component: Expandable,
        allow_extra_filters: bool = False,
    ):
        filter_dict: dict[str, list[str] | str] = {}
        filters = data.draw(
            st.lists(
                st.one_of(
                    [
                        st.sampled_from(list(component.zip_lists)),
                        st.text(string.ascii_letters, min_size=1, max_size=5),
                    ]
                )
                if allow_extra_filters
                else st.sampled_from(list(component.zip_lists)),
                max_size=len(component.zip_lists),
                unique=True,
            )
        )

        entities = {key: list(set(vals)) for key, vals in component.zip_lists.items()}

        def value_strat(filt: str):
            rand_text = st.text(sb_st.alphanum, min_size=1, max_size=10)
            return (
                st.one_of([st.sampled_from(entities[filt]), rand_text])
                if filt in entities
                else rand_text
            )

        for filt in filters:
            filter_dict[filt] = data.draw(
                st.one_of(
                    [
                        st.lists(
                            value_strat(filt),
                            unique=True,
                            min_size=1,
                            max_size=5,
                        ),
                        value_strat(filt),
                    ]
                )
            )
        return filter_dict

    def filter_iter(self, filters: dict[str, str | list[str]]):
        return {
            key: iter(val) if isinstance(val, list) else val
            for key, val in filters.items()
        }

    @given(
        component=sb_st.bids_components(max_values=4, restrict_patterns=True),
        data=st.data(),
    )
    def test_only_filter_values_in_output(
        self, component: Expandable, data: st.DataObject
    ):
        filter_dict = self.get_filter_dict(data, component)
        filtered = component.filter(**self.filter_iter(filter_dict))
        for filt in filter_dict:
            for val in filtered.zip_lists[filt]:
                assert val in filter_dict[filt]

    @given(
        component=sb_st.expandables(max_values=4, restrict_patterns=True),
        data=st.data(),
    )
    def test_zip_lists_rows_remain_of_equal_length(
        self, component: Expandable, data: st.DataObject
    ):
        filter_dict = self.get_filter_dict(data, component)
        filtered = component.filter(**self.filter_iter(filter_dict))
        lengths: set[int] = set()
        for row in filtered.zip_lists.values():
            lengths.add(len(row))
        assert len(lengths) == 1

    @given(
        component=sb_st.expandables(max_values=4, restrict_patterns=True),
        data=st.data(),
    )
    def test_all_columns_found_in_original_zip_list(
        self, component: Expandable, data: st.DataObject
    ):
        filter_dict = self.get_filter_dict(data, component, allow_extra_filters=True)
        filtered = component.filter(**self.filter_iter(filter_dict))
        cols = set(zip(*component.zip_lists.values()))
        for col in zip(*filtered.zip_lists.values()):
            assert col in cols

    @given(
        component=sb_st.expandables(max_values=4, restrict_patterns=True),
        data=st.data(),
    )
    def test_all_entities_remain_after_filtering(
        self, component: Expandable, data: st.DataObject
    ):
        filter_dict = self.get_filter_dict(data, component, allow_extra_filters=True)
        filtered = component.filter(**self.filter_iter(filter_dict))
        assert set(component.zip_lists) == set(filtered.zip_lists)

    @given(
        component=sb_st.expandables(max_values=4, restrict_patterns=True),
        data=st.data(),
    )
    def test_no_columns_that_should_be_present_are_missing(
        self,
        component: Expandable,
        data: st.DataObject,
    ):
        filter_dict = self.get_filter_dict(data, component)
        filtered = component.filter(**self.filter_iter(filter_dict))
        keys = list(component.zip_lists)
        cols = set(zip(*component.zip_lists.values()))
        result = set(zip(*filtered.zip_lists.values()))
        for col in cols:
            should_be_present = True
            for filt in filter_dict:
                if col[keys.index(filt)] not in list(
                    itx.always_iterable(filter_dict[filt])
                ):
                    assert col not in result
                    should_be_present = False
                    break
            if should_be_present:
                assert col in result


class TestFilteringBidsComponentRowWithSpec:
    def get_filter_spec(
        self,
        data: st.DataObject,
        component: BidsComponentRow,
    ):
        rand_text = st.text(sb_st.alphanum, min_size=1, max_size=10)
        value_strat = st.one_of(st.sampled_from(component.entities), rand_text)

        return data.draw(
            st.one_of(
                st.lists(
                    value_strat,
                    unique=True,
                    min_size=1,
                    max_size=5,
                ),
                value_strat,
            )
        )

    @given(
        component=sb_st.bids_component_row(max_values=4, restrict_patterns=True),
        data=st.data(),
    )
    def test_only_filter_values_in_output(
        self, component: BidsComponentRow, data: st.DataObject
    ):
        spec = self.get_filter_spec(data, component)
        filtered = component.filter(spec)
        for val in filtered.entities:
            assert val in list(itx.always_iterable(spec))

    @given(
        component=sb_st.bids_component_row(max_values=4, restrict_patterns=True),
        data=st.data(),
    )
    def test_all_columns_found_in_original_zip_list(
        self, component: BidsComponentRow, data: st.DataObject
    ):
        spec = self.get_filter_spec(data, component)
        filtered = component.filter(spec)
        cols = set(zip(*component.zip_lists.values()))
        for col in zip(*filtered.zip_lists.values()):
            assert col in cols

    @given(
        component=sb_st.bids_component_row(max_values=4, restrict_patterns=True),
        data=st.data(),
    )
    def test_all_entities_remain_after_filtering(
        self, component: BidsComponentRow, data: st.DataObject
    ):
        spec = self.get_filter_spec(data, component)
        filtered = component.filter(spec)
        assert set(component.zip_lists) == set(filtered.zip_lists)

    @given(
        component=sb_st.bids_component_row(max_values=4, restrict_patterns=True),
        data=st.data(),
    )
    def test_all_valid_values_in_spec_in_result(
        self, component: BidsComponentRow, data: st.DataObject
    ):
        spec = self.get_filter_spec(data, component)
        filtered = component.filter(itx.always_iterable(spec))
        for val in itx.always_iterable(spec):
            if val in component.entities:
                assert val in filtered.entities

    @given(
        component=sb_st.bids_component_row(max_values=4, restrict_patterns=True),
        data=st.data(),
    )
    def test_providing_both_spec_and_filters_gives_error(
        self, component: BidsComponentRow, data: st.DataObject
    ):
        spec = self.get_filter_spec(data, component)
        with pytest.raises(ValueError) as err:
            component.filter(spec, foo="bar")
        assert "__spec and filters cannot" in err.value.args[0]


class TestBidsComponentIndexing:
    def get_selectors(
        self,
        data: st.DataObject,
        dicts: BidsPartialComponent,
        use_nonexistant_keys: bool = False,
        unique: bool = False,
    ) -> tuple[str, ...]:
        if not dicts.zip_lists and not use_nonexistant_keys:
            return tuple()
        sampler = st.sampled_from(list(dicts.zip_lists))
        val_strat = (
            st.text().filter(lambda s: s not in dicts.zip_lists)
            if use_nonexistant_keys
            else sampler
        )
        return tuple(data.draw(st.lists(val_strat, unique=unique)))

    @given(dicts=sb_st.bids_partial_components(), data=st.data())
    def test_multiple_selection_returns_BidsPartialComponent(
        self, dicts: BidsPartialComponent, data: st.DataObject
    ):
        selectors = self.get_selectors(data, dicts)
        assert type(dicts[selectors]) == BidsPartialComponent

    @given(dicts=sb_st.bids_partial_components(), data=st.data())
    def test_single_selection_returns_BidsComponentRow(
        self, dicts: BidsPartialComponent, data: st.DataObject
    ):
        selectors = data.draw(st.sampled_from(list(dicts.zip_lists)))
        assert type(dicts[selectors]) == BidsComponentRow

    @given(dicts=sb_st.bids_partial_components(), data=st.data())
    def test_selected_items_in_original(
        self, dicts: BidsPartialComponent, data: st.DataObject
    ):
        selectors = self.get_selectors(data, dicts)
        selected = dicts[selectors]
        for key in selected.zip_lists:
            assert selected[key] == dicts[key]

    @given(dicts=sb_st.bids_partial_components(), data=st.data())
    def test_all_requested_items_received_and_no_others(
        self, dicts: BidsPartialComponent, data: st.DataObject
    ):
        selectors = self.get_selectors(data, dicts)
        selected = dicts[selectors]
        assert set(selectors) == set(selected.zip_lists)
        for selector in selectors:
            assert selector in selected.zip_lists

    @given(
        dicts=sb_st.bids_partial_components(),
        data=st.data(),
    )
    def test_single_missing_key_raises_error(
        self, dicts: BidsPartialComponent, data: st.DataObject
    ):
        selector = data.draw(st.text().filter(lambda s: s not in dicts.zip_lists))
        with pytest.raises(KeyError):
            dicts[selector]

    @given(dicts=sb_st.bids_partial_components(), data=st.data())
    def test_multiple_missing_key_raises_error(
        self, dicts: BidsPartialComponent, data: st.DataObject
    ):
        selectors = self.get_selectors(data, dicts, use_nonexistant_keys=True)
        assume(len(selectors))
        with pytest.raises(KeyError):
            dicts[selectors]

    @given(dicts=sb_st.bids_partial_components(), data=st.data())
    def test_order_of_selectors_is_preserved(
        self, dicts: BidsPartialComponent, data: st.DataObject
    ):
        selectors = self.get_selectors(data, dicts)
        # Using the itx method for uniqueness to avoid calculating unique values in the
        # test in the same way as the source code
        assert tuple(dicts.zip_lists[selectors]) == tuple(
            itx.unique_everseen(selectors)
        )
