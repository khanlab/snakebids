"""Helper functions and classes for tests"""

from __future__ import annotations

import functools as ft
import itertools as it
from pathlib import Path
from typing import Any, Callable, Dict, Iterable, List, TypeVar

import pytest
from hypothesis import HealthCheck, settings

from snakebids import bids
from snakebids.core.input_generation import BidsDataset, generate_inputs
from snakebids.types import InputsConfig, UserDictPy37, ZipLists
from snakebids.utils.utils import BidsEntity, MultiSelectDict

T = TypeVar("T")


def get_zip_list(
    entities: Iterable[BidsEntity | str], combinations: Iterable[tuple[str, ...]]
) -> ZipLists:
    """Return a zip list from iterables of entities and value combinations

    Parameters
    ----------
    entities : Iterable[str]
        Iterable of entities
    combinations : Iterable[tuple[str, ...]]
        Iterable of combinations (e.g. produced using itertools.product). Each tuple
        must be the same length as entities

    Returns
    -------
    dict[str, list[str]]
        zip_list representation of entity-value combinations
    """
    return MultiSelectDict(
        {
            BidsEntity(str(entity)).wildcard: list(combs)
            for entity, combs in zip(entities, zip(*combinations))
        }
    )


def setify(dic: dict[Any, list[Any]]) -> dict[Any, set[Any]]:
    """Convert a dict of *->lists into a dict of *->sets

    Parameters
    ----------
    dic : dict[Any, list[Any]]
        Dict of lists to convert

    Returns
    -------
    dict[Any, Set[Any]]
        Dict of sets
    """
    return {key: set(val) for key, val in dic.items()}


def get_bids_path(entities: Iterable[str]) -> str:
    """Get consistently ordered bids path for a group of entities

    Parameters
    ----------
    entities : Iterable[str]
        Entities to convert into a bids path

    Returns
    -------
    str
        bids path
    """

    def get_tag(entity: BidsEntity) -> tuple[str, str]:
        # For pybids, suffixes MUST be followed by extensions, but we don't yet support
        # separate indexing of extensions, so add a dummy extension any time there's a
        # suffix
        extension = ".foo" if entity == "suffix" else ""
        return entity.wildcard, f"{{{entity.wildcard}}}{extension}"

    return bids(
        root=".", **dict(get_tag(BidsEntity(entity)) for entity in sorted(entities))
    )


# pylint: disable=too-few-public-methods
class BidsListCompare(UserDictPy37[str, Dict[str, List[str]]]):
    """Dict override specifically for comparing input_lists

    When comparing, all lists are converted into sets so that order doesn't matter for
    the comparison
    """

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, dict):
            return False
        for name, lists in other.items():
            if name not in self:
                return False
            for key, val in lists.items():
                if (
                    not isinstance(val, list)
                    or key not in self[name]
                    or set(self[name][key]) != set(val)
                ):
                    return False
        return True


def debug(**overrides: Any):
    """Disable a hypothesis decorated test with specific parameter examples

    Should be used as a decorator, placed *before* @hypothesis.given(). Adds the "debug"
    mark to the test, which can be detected by other constructs (e.g. to disable fakefs
    when hypothesis is not being used)
    """

    def inner(func: Callable[..., Any]):
        test = getattr(func, "hypothesis").inner_test

        @pytest.mark.disable_fakefs(True)
        @ft.wraps(func)
        def inner_test(*args, **kwargs):
            return test(*args, **{**kwargs, **overrides})

        if not hasattr(func, "hypothesis"):
            raise TypeError(f"{func} is not decorated with hypothesis.given")
        return inner_test

    return inner


def create_dataset(root: str | Path, dataset: BidsDataset) -> None:
    """Create an empty BidsDataset on the filesystem

    Creates the directory structure and files represented by a BidsDataset. Files are
    touched: they will have no contents.
    """
    for component in dataset.values():
        entities = list(component.zip_lists.keys())
        for values in zip(*component.zip_lists.values()):
            path = Path(root, component.path.format(**dict(zip(entities, values))))
            path.parent.mkdir(parents=True, exist_ok=True)
            path.touch()


def create_snakebids_config(dataset: BidsDataset) -> InputsConfig:
    """Generate a basic snakebids config dict from a dataset"""
    all_entities = set(
        it.chain.from_iterable(comp.entities for comp in dataset.values())
    )
    return {
        comp.name: {
            "filters": {
                BidsEntity.from_tag(entity).entity: comp.entities[entity]
                if entity in comp.entities
                else False
                for entity in all_entities
            },
            "wildcards": [
                BidsEntity.from_tag(entity).entity for entity in comp.entities
            ],
        }
        for comp in dataset.values()
    }


def reindex_dataset(
    root: str, dataset: BidsDataset, use_custom_paths: bool = False
) -> BidsDataset:
    """Create BidsDataset on the filesystem and reindex

    Paths within the dataset must be absolute
    """
    create_dataset(Path("/"), dataset)
    config = create_snakebids_config(dataset)
    if use_custom_paths:
        for comp in config:
            config[comp]["custom_path"] = dataset[comp].path
    return generate_inputs(root, config, use_bids_inputs=True)


def allow_tmpdir(__callable: T) -> T:
    """Allow function_scoped fixtures in hypothesis tests

    This is primarily useful for using tmpdirs, hence, the name
    """
    return settings(
        suppress_health_check=[
            HealthCheck.function_scoped_fixture,
        ],
    )(__callable)
