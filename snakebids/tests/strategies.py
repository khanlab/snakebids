from __future__ import annotations

import copy
import itertools as it
from pathlib import Path
from string import ascii_letters, digits
from typing import Any, Optional, Type, TypeVar

import hypothesis.strategies as st
from bids.layout import Config as BidsConfig
from hypothesis import assume

from snakebids.core.datasets import BidsComponent
from snakebids.core.input_generation import BidsDataset
from snakebids.tests import helpers
from snakebids.types import InputConfig, InputsConfig, ZipList
from snakebids.utils.utils import BidsEntity, MultiSelectDict

_Ex_co = TypeVar("_Ex_co", bound=str, covariant=True)
_T = TypeVar("_T")

alphanum = ascii_letters + digits
valid_entities: tuple[str] = tuple(BidsConfig.load("bids").entities.keys())


def bids_entity() -> st.SearchStrategy[BidsEntity]:
    bidsconfig = BidsConfig.load("bids")
    # Generate inputs and bids does not properly handle 'extension', so exclude it
    return st.sampled_from(
        [
            BidsEntity(key)
            for key in bidsconfig.entities.keys()
            if key not in ["extension", "fmap", "scans"]
        ],
    )


def bids_value(pattern: str = r"[^\n\r]*") -> st.SearchStrategy[str]:
    return st.from_regex(pattern, fullmatch=True).filter(len)


def bids_entity_lists(
    min_size: int = 1, max_size: int = 5
) -> st.SearchStrategy[list[BidsEntity]]:
    return st.lists(
        bids_entity(),
        min_size=min_size,
        max_size=max_size,
        unique=True,
        # bids_paths aren't formed correctly if only datatype is provided
    ).filter(lambda v: v != ["datatype"])


@st.composite
def input_configs(draw: st.DrawFn) -> InputConfig:
    filtered_entities = draw(st.one_of(st.lists(bids_entity()), st.none()))
    filters = (
        {
            entity.entity: draw(
                st.one_of(st.booleans(), bids_value(), st.lists(bids_value()))
            )
            for entity in filtered_entities
        }
        if filtered_entities is not None
        else None
    )
    wildcard_entities = draw(st.one_of(st.lists(bids_entity()), st.none()))

    wildcards = (
        [entity.entity for entity in wildcard_entities]
        if wildcard_entities is not None
        else None
    )
    custom_path = draw(st.one_of(st.text(), st.none()))

    pybids_inputs: InputConfig = {}
    if wildcards is not None:
        pybids_inputs.update({"wildcards": wildcards})
    if filters is not None:
        pybids_inputs.update({"filters": filters})
    if custom_path is not None:
        pybids_inputs.update({"custom_path": custom_path})
    return pybids_inputs


def inputs_configs() -> st.SearchStrategy[InputsConfig]:
    return st.dictionaries(st.text(min_size=1), input_configs())


@st.composite
def zip_lists(  # noqa: PLR0913
    draw: st.DrawFn,
    min_entities: int = 1,
    max_entities: int = 5,
    min_values: int = 1,
    max_values: int = 3,
    entities: Optional[list[BidsEntity]] = None,
    restrict_patterns: bool = False,
) -> ZipList:
    # Generate multiple entity sets for different "file types"

    if entities is None:
        entities = draw(bids_entity_lists(min_size=min_entities, max_size=max_entities))

    values = {
        entity: draw(
            st.lists(
                bids_value(entity.match if restrict_patterns else ".*"),
                min_size=min_values,
                max_size=max_values,
                unique=True,
            )
        )
        for entity in entities
    }

    combinations = list(it.product(*values.values()))
    used_combinations = draw(
        st.lists(
            st.sampled_from(combinations),
            min_size=1,
            max_size=len(combinations),
            unique=True,
        )
    )
    return helpers.get_zip_list(values, used_combinations)


@st.composite
def bids_components(  # noqa: PLR0913
    draw: st.DrawFn,
    min_entities: int = 1,
    max_entities: int = 5,
    min_values: int = 1,
    max_values: int = 3,
    entities: Optional[list[BidsEntity]] = None,
    root: Optional[Path] = None,
    name: str | None = None,
    restrict_patterns: bool = False,
) -> BidsComponent:
    zip_list = draw(
        zip_lists(
            min_entities=min_entities,
            max_entities=max_entities,
            min_values=min_values,
            max_values=max_values,
            entities=entities,
            restrict_patterns=restrict_patterns,
        )
    )

    path = (root or Path()) / helpers.get_bids_path(zip_list)

    return BidsComponent(
        name=name or draw(bids_value()),
        path=str(path),
        zip_lists=zip_list,
    )


@st.composite
def bids_input_lists(
    draw: st.DrawFn,
    min_size: int = 1,
    max_size: int = 5,
    entities: Optional[list[BidsEntity]] = None,
) -> dict[str, list[str]]:
    # Generate multiple entity sets for different "file types"
    if entities is None:
        entities = draw(bids_entity_lists(min_size))

    return {
        entity.wildcard: draw(
            st.lists(bids_value(), min_size=min_size, max_size=max_size, unique=True)
        )
        for entity in entities
    }


def everything_except(*excluded_types: Type[Any]) -> st.SearchStrategy[Any]:
    return (
        st.from_type(type)
        .flatmap(st.from_type)
        .filter(lambda s: not isinstance(s, excluded_types))
    )


@st.composite
def datasets(
    draw: st.DrawFn,
    root: Optional[Path] = None,
) -> BidsDataset:
    ent1 = draw(bids_entity_lists(min_size=2, max_size=3))
    ent2 = copy.copy(ent1)
    ent2.pop()
    # BUG: need better controlling if 'datatype' is the only arg
    assume(ent2 != ["datatype"])
    comp1 = draw(bids_components(entities=ent1, restrict_patterns=True, root=root))
    comp2 = draw(bids_components(entities=ent2, restrict_patterns=True, root=root))
    assume(comp1.input_name != comp2.input_name)
    return BidsDataset.from_iterable([comp1, comp2])


@st.composite
def datasets_one_comp(
    draw: st.DrawFn,
    root: Optional[Path] = None,
    names: st.SearchStrategy[str] | None = None,
) -> BidsDataset:
    ent1 = draw(bids_entity_lists(min_size=2, max_size=3))
    comp1 = draw(
        bids_components(
            entities=ent1,
            restrict_patterns=True,
            root=root,
            name=draw(names) if names is not None else None,
        )
    )
    return BidsDataset.from_iterable([comp1])


@st.composite
def multiselect_dicts(
    draw: st.DrawFn,
    keys: st.SearchStrategy[_Ex_co],
    values: st.SearchStrategy[_T],
    *,
    min_size: int = 0,
    max_size: Optional[int] = None,
) -> MultiSelectDict[_Ex_co, _T]:
    return MultiSelectDict(
        draw(
            st.dictionaries(
                keys,
                values,
                min_size=min_size,
                max_size=max_size,
            )
        )
    )


def everything() -> st.SearchStrategy[Any]:
    return st.from_type(type).flatmap(st.from_type)
