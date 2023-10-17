"""Helper functions and classes for tests"""
from __future__ import annotations

import functools as ft
import itertools as it
import re
import subprocess as sp
from datetime import timedelta
from pathlib import Path
from typing import (
    Any,
    Callable,
    Dict,
    Iterable,
    List,
    Mapping,
    Protocol,
    Sequence,
    TypeVar,
)

import more_itertools as itx
import pytest
from hypothesis import HealthCheck, example, settings
from typing_extensions import ParamSpec

from snakebids import bids
from snakebids.core.datasets import BidsDataset
from snakebids.core.input_generation import generate_inputs
from snakebids.types import InputsConfig, UserDictPy38, ZipList, ZipListLike
from snakebids.utils.utils import BidsEntity, MultiSelectDict

_T = TypeVar("_T")
_T_contra = TypeVar("_T_contra", contravariant=True)
_P = ParamSpec("_P")


def get_zip_list(
    entities: Iterable[BidsEntity | str], combinations: Iterable[tuple[str, ...]]
) -> ZipList:
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

    def strlist() -> list[str]:
        return []

    lists: Iterable[Sequence[str]] = list(zip(*combinations)) or itx.repeatfunc(strlist)
    return MultiSelectDict(
        {
            BidsEntity(str(entity)).wildcard: list(combs)
            for entity, combs in zip(entities, lists)
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


def get_bids_path(entities: Iterable[str | BidsEntity], **extras: str) -> str:
    """Get consistently ordered bids path for a group of entities

    Parameters
    ----------
    entities : Iterable[str]
        Entities to convert into a bids path
    extras
        Extra entity-value pairs to include in the path

    Returns
    -------
    str
        bids path
    """

    def get_tag(entity: BidsEntity) -> tuple[str, str]:
        return entity.wildcard, f"{{{entity.wildcard}}}"

    return bids(
        **dict(get_tag(BidsEntity(entity)) for entity in sorted(entities)),
        **{BidsEntity(entity).wildcard: value for entity, value in extras.items()},
    )


class BidsListCompare(UserDictPy38[str, Dict[str, List[str]]]):
    """Dict override specifically for comparing input_lists

    When comparing, all lists are converted into sets so that order doesn't matter for
    the comparison
    """

    def __eq__(self, other: object) -> bool:
        # Getting pyright to like this is more work than it's worth
        if not isinstance(other, dict):
            return False
        for name, lists in other.items():  # type: ignore
            if name not in self:
                return False
            for key, val in lists.items():  # type: ignore
                if (
                    not isinstance(val, list)
                    or key not in self[name]  # type: ignore
                    or set(self[name][key]) != set(val)  # type: ignore
                ):
                    return False
        return True


def debug(**overrides: Any):
    """Disable a hypothesis decorated test with specific parameter examples

    Should be used as a decorator, placed *before* @hypothesis.given(). Adds the "debug"
    mark to the test, which can be detected by other constructs (e.g. to disable fakefs
    when hypothesis is not being used)
    """

    def inner(func: Callable[_P, _T]) -> Callable[_P, _T]:
        if not hasattr(func, "hypothesis"):
            msg = f"{func} is not decorated with hypothesis.given"
            raise TypeError(msg)

        test = func.hypothesis.inner_test  # type: ignore

        @pytest.mark.disable_fakefs(True)
        @ft.wraps(func)
        def inner_test(*args: _P.args, **kwargs: _P.kwargs):
            return test(*args, **{**kwargs, **overrides})

        return inner_test

    return inner


def mock_data(*draws: Any):
    """Utility function for mocking the hypothesis data strategy

    Intended for combination with debug. Takes a list of values corresponding the draws
    taking from the data object. In the debug run, calls to data will return the given
    values in the order specified
    """

    # pylint: disable=missing-class-docstring, too-few-public-methods
    class MockData:
        _draws = iter(draws)

        # pylint: disable=unused-argument
        def draw(self, strategy: Any, label: Any = None):
            return next(self._draws)

    return MockData()


def create_dataset(root: str | Path, dataset: BidsDataset) -> None:
    """Create an empty BidsDataset on the filesystem

    Creates the directory structure and files represented by a BidsDataset. Files are
    touched: they will have no contents.
    """
    for component in dataset.values():
        for path in map(Path(root).joinpath, component.expand()):
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
                BidsEntity.from_tag(entity).entity: list(comp.entities[entity])
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
    root: str,
    dataset: BidsDataset,
    use_custom_paths: bool = False,
    participant_label: str | Sequence[str] | None = None,
    exclude_participant_label: str | Sequence[str] | None = None,
) -> BidsDataset:
    """Create BidsDataset on the filesystem and reindex

    Paths within the dataset must be absolute
    """
    create_dataset(Path("/"), dataset)
    config = create_snakebids_config(dataset)
    if participant_label is not None or exclude_participant_label is not None:
        for comp in config.values():
            if "subject" in comp.get("filters", {}):
                del comp["filters"]["subject"]  # type: ignore

    if use_custom_paths:
        for comp in config:
            config[comp]["custom_path"] = dataset[comp].path
    return generate_inputs(
        root,
        config,
        participant_label=participant_label,
        exclude_participant_label=exclude_participant_label,
    )


def allow_function_scoped(func: _T, /) -> _T:
    """Allow function_scoped fixtures in hypothesis tests

    This is primarily useful for using tmpdirs, hence, the name
    """
    return settings(
        suppress_health_check=[
            HealthCheck.function_scoped_fixture,
        ],
    )(func)


def deadline(time: int | float | timedelta | None) -> Callable[[_T], _T]:
    """Change hypothesis deadline

    Numbers refer to time in milliseconds. Set to None to disable entirely
    """

    def inner(func: _T, /) -> _T:
        return settings(deadline=time)(func)

    return inner


def expand_zip_list(
    zip_list: ZipListLike, new_values: Mapping[str, Sequence[str]]
) -> ZipListLike:
    """Expand a zip list with new values via product

    Each zip list column will be combined with all possible combinations of the
    entity-value pairs in ``new_values``. It is thus equivalent to calling
    :meth:`BidsComponent.expand() <snakebids.BidsComponent.expand>` with ``new_values``
    as extra args.
    """
    zip_cols = list(zip(*zip_list.values()))
    new_cols = list(
        map(
            it.chain.from_iterable,
            it.product(zip_cols, it.product(*new_values.values())),
        )
    )
    z = zip(*new_cols)
    return dict(zip(it.chain(zip_list.keys(), new_values.keys()), z))


def needs_docker(container: str):
    def decorator(func: Callable[_P, _T]) -> Callable[_P, _T]:
        @pytest.mark.docker()
        @ft.wraps(func)
        def wrapper(*args: _P.args, **kwargs: _P.kwargs):
            try:
                sp.run(["docker"], check=True)
            except (sp.CalledProcessError, FileNotFoundError):
                pytest.fail(
                    "docker is not available on this machine. To skip docker tests, "
                    "use '-m \"not docker\"'"
                )
            try:
                sp.run(["docker", "image", "inspect", container], check=True)
            except sp.CalledProcessError as err:
                if not (
                    match := re.match(r"snakebids/([\w\-])+\:[a-zA-Z0-9\-]+", container)
                ):
                    pytest.fail(
                        f"Unrecognized docker image: {container}. Should be "
                        "'snakebids/{container_id}:{version}"
                    )
                container_id = match[1]
                pytest.fail(
                    f"{container} is not built on this machine. To build container, "
                    f"run `poetry run poe build-container {container_id}`. To skip "
                    f"docker tests, use '-m \"not docker\"'. (got error {err})"
                )
            return func(*args, **kwargs)

        return wrapper

    return decorator


def entity_to_wildcard(entities: str | Iterable[str], /):
    """Turn entity strings into wildcard dicts as {"entity": "{entity}"}"""
    return {entity: f"{{{entity}}}" for entity in itx.always_iterable(entities)}


def identity(obj: _T) -> _T:
    return obj


def example_if(condition: bool, *args: Any, **kwargs: Any):
    def inner(func: Callable[_P, _T]) -> Callable[_P, _T]:
        return example(*args, **kwargs)(func)

    if condition:
        return inner

    return identity


class Benchmark(Protocol):
    def __call__(
        self, func: Callable[_P, _T], *args: _P.args, **kwargs: _P.kwargs
    ) -> _T:
        ...


"""Comparison Dunders copied from typeshed"""


class SupportsDunderLT(Protocol[_T_contra]):
    def __lt__(self, __other: _T_contra) -> bool:
        ...


def is_strictly_increasing(items: Iterable[SupportsDunderLT[Any]]) -> bool:
    # itx.pairwise properly aliases it.pairwise on py310+
    return all(i < j for i, j in itx.pairwise(items))
