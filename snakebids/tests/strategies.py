import itertools as it
from typing import Any, List, Optional, Type

import hypothesis.strategies as st
from bids.layout import Config as BidsConfig

from snakebids import BidsComponent
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
                bids_value(entity.pattern if restrict_patterns else ".*"),
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

    path = helpers.get_bids_path(zip_lists)

    return BidsComponent(
        name=draw(bids_value()),
        path=path,
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
