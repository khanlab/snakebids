from __future__ import annotations

import os
from collections.abc import Iterable
from pathlib import Path

import pytest
from hypothesis import given
from hypothesis import strategies as st
from pathvalidate import Platform, is_valid_filename

from snakebids.paths.presets import bids as bids_old
from snakebids.paths.presets import bids_b as bids
from snakebids.paths.specs import v0_0_0
from snakebids.tests import strategies as sb_st
from snakebids.tests.helpers import Benchmark, is_strictly_increasing
from snakebids.utils.utils import BidsEntity

V0_0_0 = v0_0_0()


def _get_entity_tags(entities: Iterable[str]):
    return [BidsEntity.normalize(entity).tag for entity in entities]


def _values():
    return st.text(min_size=1).filter(
        lambda s: is_valid_filename(s, Platform.LINUX) and "_" not in s and "-" not in s
    )


def _bids_args(schema_only: bool = True):
    has_tag = [e["entity"] for e in V0_0_0 if e.get("tag")]
    return (
        st.dictionaries(
            keys=st.one_of(
                sb_st.bids_entity(whitelist_entities=[e["entity"] for e in V0_0_0]),
                sb_st.bids_entity(
                    whitelist_entities=["datatype", "suffix", "extension"]
                )
                if not schema_only
                else sb_st.nothing(),
            ),
            # The boolean here is to decide whether to use the entity or the tag in the
            # BidsEntity generated above
            values=st.tuples(
                st.booleans(),
                _values(),
            ),
            min_size=1,
        )
        .map(
            # we only want to use the full entity name if the spec in use understands
            # the full entity name
            lambda d: {
                key.entity if val[0] and str(key) in has_tag else key.tag: val[1]
                for key, val in d.items()
            }
        )
        .filter(lambda s: set(s) - {"datatype", "extension"})
    )


@given(_bids_args())
def test_number_of_underscore_corresponds_to_number_entities(entities: dict[str, str]):
    assert bids(**entities).count("_") == len(entities) - 1


@given(_bids_args())
def test_number_of_dashes_corresponds_to_number_entities(entities: dict[str, str]):
    assert Path(bids(**entities)).name.count("-") == len(entities)


@given(entities=_bids_args(schema_only=False), prefix=_values())
def test_beginning_of_name_always_prefix(entities: dict[str, str], prefix: str):
    assert Path(bids(prefix=prefix, **entities)).name.startswith(prefix)


@given(entities=_bids_args(schema_only=False), root=_values())
def test_beginning_of_path_always_root(entities: dict[str, str], root: str):
    path = bids(root=root, **entities)
    assert path.startswith(root)
    assert path[len(root)] == os.path.sep


@given(entities=_bids_args(), datatype=_values())
def test_one_path_element_always_datatype(entities: dict[str, str], datatype: str):
    assert datatype in {
        p.name for p in Path(bids(datatype=datatype, **entities)).parents
    }


@given(entities=_bids_args())
def test_entities_all_in_path_as_tags(entities: dict[str, str]):
    tags = _get_entity_tags(entities)
    path = bids(**entities)
    for tag in tags:
        assert f"{tag}-" in path


@given(entities=_bids_args())
def test_full_entity_names_not_in_path(entities: dict[str, str]):
    non_tags = [
        normalized
        for entity in entities
        if (normalized := BidsEntity.normalize(entity).tag) != entity
    ]
    path = bids(**entities)
    for entity in non_tags:
        assert f"{entity}-" not in path


@given(entities=_bids_args())
def test_entities_found_in_name_in_correct_order(entities: dict[str, str]):
    tags = _get_entity_tags(entities)
    order: list[str] = []
    for entity in V0_0_0:
        if (tag := entity.get("tag", entity["entity"])) in tags:
            order.append(tag)
    path = bids(**entities)
    assert set(tags) == set(order)
    assert is_strictly_increasing(map(path.index, order))


def test_bids_subj():
    assert bids(root="bids", subject="001", suffix="T1w.nii.gz") == (
        "bids/sub-001/sub-001_T1w.nii.gz"
    )
    assert bids(root=Path("bids").resolve(), subject="001", suffix="T1w.nii.gz") == (
        str(Path.cwd() / "bids/sub-001/sub-001_T1w.nii.gz")
    )


def test_bids_with_no_args_gives_empty_path():
    assert not bids()


@pytest.mark.parametrize(
    ("args"),
    [
        {"datatype": "foo"},
        {"prefix": "foo"},
        {"datatype": "foo", "prefix": "foo"},
    ],
)
def test_missing_essential_entities_gives_error(args: dict[str, str]):
    with pytest.raises(ValueError):
        bids(**args)


def test_if_original_is_faster(benchmark: Benchmark):
    benchmark(
        bids_old,
        root="foo/bar",
        subject="001",
        session="32",
        run="first",
        space="mine",
        foo="bar",
        england="britain",
        rome="duckling",
    )


def test_if_new_is_faster(benchmark: Benchmark):
    benchmark(
        bids,
        root="foo/bar",
        subject="001",
        session="32",
        run="first",
        space="mine",
        foo="bar",
        england="britain",
        rome="duckling",
    )
