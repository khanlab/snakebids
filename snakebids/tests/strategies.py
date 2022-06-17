import itertools as it
from typing import Any, Dict, List, Type, cast

import hypothesis.strategies as st
from bids.layout import Config as BidsConfig

from snakebids import BidsComponent, bids
from snakebids.tests.helpers import get_zip_list


def bids_entity():
    bidsconfig = BidsConfig.load("bids")
    return st.sampled_from(cast(str, list(bidsconfig.entities.keys())))


def bids_value():
    bids_alphabet = st.characters(whitelist_categories=["Ll", "Lu", "Nd"])
    return st.text(bids_alphabet, min_size=1, max_size=10)


def bids_entity_list(min_size: int = 1):
    return st.lists(
        bids_entity(),
        min_size=min_size,
        max_size=5,
        unique=True,
    )


def get_bids_path(zip_lists: Dict[str, List[str]]):
    return bids(root=".", **{entity: f"{{{entity}}}" for entity in zip_lists})


@st.composite
def bids_input(draw: st.DrawFn):
    # Generate multiple entity sets for different "file types"

    entities = draw(bids_entity_list())

    values = {
        key: draw(st.lists(bids_value(), min_size=1, unique=True)) for key in entities
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
    zip_lists = get_zip_list(values, used_combinations)

    path = get_bids_path(zip_lists)

    return BidsComponent(
        input_name=draw(bids_value()),
        input_path=path,
        input_zip_lists=zip_lists,
    )


@st.composite
def bids_input_lists(draw: st.DrawFn, min_size: int = 1, max_size: int = 5):
    # Generate multiple entity sets for different "file types"
    entities = draw(bids_entity_list(min_size))

    return {
        key: draw(
            st.lists(bids_value(), min_size=min_size, max_size=max_size, unique=True)
        )
        for key in entities
    }


def everything_except(*excluded_types: Type[Any]):
    return (
        st.from_type(type)
        .flatmap(st.from_type)
        .filter(lambda s: not isinstance(s, excluded_types))
    )
