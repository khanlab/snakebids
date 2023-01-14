import copy
import itertools as it
from pathlib import Path
from typing import Any, List, Optional, Type

import hypothesis.strategies as st
from bids.layout import Config as BidsConfig
from hypothesis import assume

from snakebids.core.datasets import BidsComponent
from snakebids.core.input_generation import BidsDataset
from snakebids.tests import helpers
from snakebids.utils.utils import BidsEntity


def bids_entity():
    bidsconfig = BidsConfig.load("bids")
    # Generate inputs and bids does not properly handle 'extension', so exclude it
    return st.sampled_from(
        [
            BidsEntity(key)
            for key in bidsconfig.entities.keys()
            if key not in ["extension", "fmap", "scans"]
        ],
    )


def bids_value(pattern: str = ".*"):
    return st.from_regex(pattern, fullmatch=True).filter(len)


def bids_entity_lists(min_size: int = 1, max_size: int = 5):
    return st.lists(
        bids_entity(),
        min_size=min_size,
        max_size=max_size,
        unique=True,
        # bids_paths aren't formed correctly if only datatype is provided
    ).filter(lambda v: v != ["datatype"])


@st.composite
def input_zip_lists(
    draw: st.DrawFn,
    min_size: int = 1,
    max_size: int = 5,
    entities: Optional[List[BidsEntity]] = None,
    restrict_patterns: bool = False,
):
    # Generate multiple entity sets for different "file types"

    if entities is None:
        entities = draw(bids_entity_lists(min_size=min_size, max_size=max_size))

    # TODO: min_size and max_size shouldn't be hard-coded, but we need a good way of
    # doing this.
    values = {
        entity: draw(
            st.lists(
                bids_value(entity.match if restrict_patterns else ".*"),
                min_size=1,
                max_size=3,
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
def bids_components(
    draw: st.DrawFn,
    min_size: int = 1,
    max_size: int = 5,
    entities: Optional[List[BidsEntity]] = None,
    root: Optional[Path] = None,
    restricted_chars: bool = False,
):
    zip_lists = draw(
        input_zip_lists(
            min_size=min_size,
            max_size=max_size,
            entities=entities,
            restrict_patterns=restricted_chars,
        )
    )

    path = (root or Path()) / helpers.get_bids_path(zip_lists)

    return BidsComponent(
        name=draw(bids_value()),
        path=str(path),
        zip_lists=zip_lists,
    )


@st.composite
def bids_input_lists(
    draw: st.DrawFn,
    min_size: int = 1,
    max_size: int = 5,
    entities: Optional[List[BidsEntity]] = None,
):
    # Generate multiple entity sets for different "file types"
    if entities is None:
        entities = draw(bids_entity_lists(min_size))

    return {
        entity.wildcard: draw(
            st.lists(bids_value(), min_size=min_size, max_size=max_size, unique=True)
        )
        for entity in entities
    }


def everything_except(*excluded_types: Type[Any]):
    return (
        st.from_type(type)
        .flatmap(st.from_type)
        .filter(lambda s: not isinstance(s, excluded_types))
    )


@st.composite
def datasets(draw: st.DrawFn, root: Optional[Path] = None):
    ent1 = draw(bids_entity_lists(min_size=2, max_size=3))
    assume("datatype" not in ent1)
    # Currently, space and ce cannot coexist because ce is a substr of space (see
    # snakebids.core.input_generation:_parse_bids_path)
    assume(not ("space" in ent1 and "ceagent" in ent1))
    ent2 = copy.copy(ent1)
    ent2.pop()
    # BUG: snakebids currently doesn't properly parse paths with just suffix
    assume(ent2 != ["suffix"])
    comp1 = draw(bids_components(entities=ent1, restricted_chars=True, root=root))
    comp2 = draw(bids_components(entities=ent2, restricted_chars=True, root=root))
    assume(comp1.input_name != comp2.input_name)
    return BidsDataset.from_iterable([comp1, comp2])


@st.composite
def datasets_one_comp(draw: st.DrawFn, root: Optional[Path] = None):
    ent1 = draw(bids_entity_lists(min_size=2, max_size=3))
    assume("datatype" not in ent1)
    # Currently, space and ce cannot coexist because ce is a substr of space (see
    # snakebids.core.input_generation:_parse_bids_path)
    assume(not ("space" in ent1 and "ceagent" in ent1))
    comp1 = draw(bids_components(entities=ent1, restricted_chars=True, root=root))
    return BidsDataset.from_iterable([comp1])
