from __future__ import annotations

import copy
import itertools as it
import sys
from os import PathLike
from pathlib import Path
from string import ascii_letters, digits
from typing import Any, Container, Hashable, Iterable, Optional, Type, TypeVar

import hypothesis.strategies as st
from bids.layout import Config as BidsConfig
from hypothesis import assume

from snakebids.core.datasets import (
    BidsComponent,
    BidsComponentRow,
    BidsDataset,
    BidsPartialComponent,
)
from snakebids.tests import helpers
from snakebids.types import Expandable, InputConfig, InputsConfig, ZipList
from snakebids.utils.utils import BidsEntity, MultiSelectDict

_Ex_co = TypeVar("_Ex_co", bound=str, covariant=True)
_T = TypeVar("_T")

alphanum = ascii_letters + digits
valid_entities: tuple[str] = tuple(BidsConfig.load("bids").entities.keys())


def bids_entity(
    *,
    blacklist_entities: Optional[Container[BidsEntity | str]] = None,
    whitelist_entities: Optional[Container[BidsEntity | str]] = None,
    path_safe: bool = False,
) -> st.SearchStrategy[BidsEntity]:
    blacklist = (
        helpers.ContainerBag(
            blacklist_entities if blacklist_entities is not None else set(),
            {"datatype", "suffix", "extension"},
        )
        if path_safe
        else blacklist_entities or set()
    )
    return st.sampled_from(
        [
            BidsEntity(key)
            for key in valid_entities
            if key not in ["fmap", "scans"]
            and key not in (blacklist)
            and (not whitelist_entities or key in whitelist_entities)
        ],
    )


def bids_value(pattern: str = r"[^\n\r]*") -> st.SearchStrategy[str]:
    # Need the isdigit == isdecimal check to work around a pybids bug
    return (
        st.from_regex(pattern, fullmatch=True)
        .filter(len)
        .filter(lambda s: all(c.isdigit() == c.isdecimal() for c in s))
        # remove values leading with multiple dots ".." because pybids doesn't handle
        # them correctly when applied to the "extension" entity
        .filter(lambda s: s[:2] != "..")
    )


def _filter_invalid_entity_lists(entities: Container[BidsEntity | str]):
    """Entity lists may not consist of just a datatype or extension.

    If suffix is in the entity list, so must extension
    """
    # The versions of pybids available from py37 have bugs in the indexing of suffix
    # and extension, so just exclude these entries from our tests
    if sys.version_info < (3, 8) and ("extension" in entities or "suffix" in entities):
        return False
    return all(
        [
            # If suffix is in the path, extension must be too
            ("suffix" not in entities or "extension" in entities),
            # Cannot have paths with just datatype, just extension, or just datatype and
            # extension
            entities not in [["datatype"], ["extension"], ["datatype", "extension"]],
        ]
    )


@st.composite
def bids_path(  # noqa: PLR0913
    draw: st.DrawFn,
    *,
    root: PathLike[str] | str | None = None,
    entities: Iterable[BidsEntity | str] | None = None,
    blacklist_entities: Optional[Container[BidsEntity | str]] = None,
    whitelist_entities: Optional[Container[BidsEntity | str]] = None,
    extra_entities: bool = True,
) -> Path:
    entities = (
        draw(
            bids_entity_lists(
                whitelist_entities=whitelist_entities,
                blacklist_entities=blacklist_entities,
            )
        )
        if entities is None
        else list(entities)
    )

    extras = (
        {
            k: v[0].replace("{", "{{").replace("}", "}}")
            for k, v in draw(
                zip_lists(
                    max_values=1,
                    max_entities=2,
                    blacklist_entities=helpers.ContainerBag(
                        blacklist_entities if blacklist_entities is not None else [],
                        [BidsEntity.normalize(e) for e in entities],
                    ),
                    restrict_patterns=True,
                )
            ).items()
        }
        if extra_entities
        else {}
    )

    return Path(root or ".") / helpers.get_bids_path(entities, **extras)


def bids_entity_lists(
    *,
    min_size: int = 1,
    max_size: int = 5,
    blacklist_entities: Optional[Container[BidsEntity | str]] = None,
    whitelist_entities: Optional[Container[BidsEntity | str]] = None,
) -> st.SearchStrategy[list[BidsEntity]]:
    return st.lists(
        bids_entity(
            whitelist_entities=whitelist_entities,
            blacklist_entities=blacklist_entities,
        ),
        min_size=min_size,
        max_size=max_size,
        unique=True,
    ).filter(_filter_invalid_entity_lists)


@st.composite
def input_configs(draw: st.DrawFn) -> InputConfig:
    filters = draw(
        st.one_of(
            st.dictionaries(
                bids_entity().map(str),
                st.one_of(st.booleans(), bids_value(), st.lists(bids_value())),
            ),
            st.none(),
        )
    )

    wildcards = draw(st.one_of(st.lists(bids_entity().map(str)), st.none()))
    custom_path = draw(
        st.one_of(
            st.text(
                alphabet=st.characters(
                    blacklist_categories=("Cs",), blacklist_characters=("\x00",)
                )
            ),
            st.none(),
        )
    )

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
    *,
    min_entities: int = 1,
    max_entities: int = 5,
    min_values: int = 1,
    max_values: int = 3,
    entities: Optional[list[BidsEntity]] = None,
    blacklist_entities: Optional[Container[BidsEntity | str]] = None,
    whitelist_entities: Optional[Container[BidsEntity | str]] = None,
    restrict_patterns: bool = False,
    unique: bool = False,
    cull: bool = True,
) -> ZipList:
    # Generate multiple entity sets for different "file types"

    if entities is None:
        entities = draw(
            bids_entity_lists(
                min_size=min_entities,
                max_size=max_entities,
                blacklist_entities=blacklist_entities,
                whitelist_entities=whitelist_entities,
            )
        )

    values = {
        entity: draw(
            st.lists(
                bids_value(entity.match if restrict_patterns else ".*"),
                min_size=min_values,
                max_size=max_values,
                unique=(cull or unique),
            )
        )
        for entity in entities
    }

    combinations = list(it.product(*values.values()))
    used_combinations = (
        draw(
            st.lists(
                st.sampled_from(combinations),
                min_size=1,
                max_size=len(combinations),
                unique=unique,
            )
        )
        if cull
        else combinations
    )
    return helpers.get_zip_list(values, used_combinations)


@st.composite
def bids_component_row(  # noqa: PLR0913
    draw: st.DrawFn,
    *,
    min_values: int = 1,
    max_values: int = 5,
    entity: Optional[BidsEntity] = None,
    blacklist_entities: Optional[Container[BidsEntity | str]] = None,
    whitelist_entities: Optional[Container[BidsEntity | str]] = None,
    restrict_patterns: bool = False,
    unique: bool = False,
    path_safe: bool = False,
) -> BidsComponentRow:
    entity = entity or draw(
        bids_entity(
            blacklist_entities=blacklist_entities,
            whitelist_entities=whitelist_entities,
            path_safe=path_safe,
        )
    )
    values = draw(
        st.lists(
            bids_value(entity.match if restrict_patterns else ".*"),
            min_size=min_values,
            max_size=max_values,
            unique=unique,
        )
    )
    return BidsComponentRow(values, entity=entity.entity)


@st.composite
def bids_partial_components(  # noqa: PLR0913
    draw: st.DrawFn,
    *,
    min_entities: int = 1,
    max_entities: int = 5,
    min_values: int = 1,
    max_values: int = 3,
    entities: Optional[list[BidsEntity]] = None,
    blacklist_entities: Optional[Container[BidsEntity | str]] = None,
    whitelist_entities: Optional[Container[BidsEntity | str]] = None,
    restrict_patterns: bool = False,
    unique: bool = False,
    cull: bool = True,
    _allow_subclasses: bool = True,
) -> BidsPartialComponent:
    if _allow_subclasses:
        return draw(
            st.one_of(
                bids_components(
                    min_entities=min_entities,
                    max_entities=max_entities,
                    min_values=min_values,
                    max_values=max_values,
                    entities=entities,
                    blacklist_entities=blacklist_entities,
                    whitelist_entities=whitelist_entities,
                    restrict_patterns=restrict_patterns,
                    cull=cull,
                    unique=unique,
                ),
                bids_partial_components(
                    min_entities=min_entities,
                    max_entities=max_entities,
                    min_values=min_values,
                    max_values=max_values,
                    entities=entities,
                    blacklist_entities=blacklist_entities,
                    whitelist_entities=whitelist_entities,
                    restrict_patterns=restrict_patterns,
                    cull=cull,
                    unique=unique,
                    _allow_subclasses=False,
                ),
            )
        )
    zip_list = draw(
        zip_lists(
            min_entities=min_entities,
            max_entities=max_entities,
            min_values=min_values,
            max_values=max_values,
            entities=entities,
            blacklist_entities=blacklist_entities,
            whitelist_entities=whitelist_entities,
            restrict_patterns=restrict_patterns,
            cull=cull,
            unique=unique,
        )
    )

    return BidsPartialComponent(zip_lists=zip_list)


@st.composite
def bids_components(  # noqa: PLR0913
    draw: st.DrawFn,
    *,
    min_entities: int = 1,
    max_entities: int = 5,
    min_values: int = 1,
    max_values: int = 3,
    entities: Optional[list[BidsEntity]] = None,
    blacklist_entities: Optional[Container[BidsEntity | str]] = None,
    whitelist_entities: Optional[Container[BidsEntity | str]] = None,
    root: Optional[Path] = None,
    name: str | None = None,
    restrict_patterns: bool = False,
    extra_entities: bool = True,
    blacklist_extra_entities: Optional[Container[BidsEntity | str]] = None,
    whitelist_extra_entities: Optional[Container[BidsEntity | str]] = None,
    unique: bool = False,
    cull: bool = True,
) -> BidsComponent:
    partial = draw(
        bids_partial_components(
            min_entities=min_entities,
            max_entities=max_entities,
            min_values=min_values,
            max_values=max_values,
            entities=entities,
            blacklist_entities=blacklist_entities,
            whitelist_entities=whitelist_entities,
            restrict_patterns=restrict_patterns,
            cull=cull,
            unique=unique,
            _allow_subclasses=False,
        )
    )
    path = draw(
        bids_path(
            root=root,
            entities=partial.zip_lists,
            extra_entities=extra_entities,
            blacklist_entities=blacklist_extra_entities,
            whitelist_entities=whitelist_extra_entities,
        )
    )
    return BidsComponent(
        name=name or draw(bids_value()),
        path=str(path),
        zip_lists=partial.zip_lists,
    )


@st.composite
def expandables(  # noqa: PLR0913
    draw: st.DrawFn,
    *,
    max_entities: int = 5,
    min_values: int = 1,
    max_values: int = 3,
    entities: Optional[list[BidsEntity]] = None,
    blacklist_entities: Optional[Container[BidsEntity | str]] = None,
    whitelist_entities: Optional[Container[BidsEntity | str]] = None,
    restrict_patterns: bool = False,
    unique: bool = False,
    cull: bool = True,
    path_safe: bool = False,
) -> Expandable:
    def get_entity(_entities: Optional[list[BidsEntity]]) -> BidsEntity | None:
        if not _entities:
            return None
        return draw(st.sampled_from(_entities))

    return draw(
        st.one_of(
            bids_partial_components(
                min_entities=1,
                max_entities=max_entities,
                min_values=min_values,
                max_values=max_values,
                entities=entities,
                blacklist_entities=blacklist_entities,
                whitelist_entities=whitelist_entities,
                restrict_patterns=restrict_patterns,
                cull=cull,
                unique=unique,
            ),
            bids_component_row(
                min_values=min_values,
                max_values=max_values,
                entity=get_entity(entities),
                blacklist_entities=blacklist_entities,
                whitelist_entities=whitelist_entities,
                restrict_patterns=restrict_patterns,
                unique=unique,
                path_safe=path_safe,
            ),
        )
    )


@st.composite
def bids_input_lists(
    draw: st.DrawFn,
    *,
    min_size: int = 1,
    max_size: int = 5,
    entities: Optional[list[BidsEntity]] = None,
) -> dict[str, list[str]]:
    # Generate multiple entity sets for different "file types"
    if entities is None:
        entities = draw(bids_entity_lists(min_size=min_size))

    return {
        entity.wildcard: draw(
            st.lists(bids_value(), min_size=min_size, max_size=max_size, unique=True)
        )
        for entity in entities
    }


@st.composite
def datasets(
    draw: st.DrawFn,
    *,
    root: Optional[Path] = None,
    unique: bool = False,
    cull: bool = True,
) -> BidsDataset:
    ent1 = draw(bids_entity_lists(min_size=2, max_size=3))
    ent2 = copy.copy(ent1)
    extra = ent2.pop()
    assume(_filter_invalid_entity_lists(ent2))
    comp1 = draw(
        bids_components(
            entities=ent1, restrict_patterns=True, root=root, unique=unique, cull=cull
        )
    )
    comp2 = draw(
        bids_components(
            entities=ent2,
            restrict_patterns=True,
            root=root,
            blacklist_extra_entities=[extra],
            unique=unique,
            cull=cull,
        )
    )
    assume(comp1.input_name != comp2.input_name)
    return BidsDataset.from_iterable([comp1, comp2])


@st.composite
def datasets_one_comp(  # noqa: PLR0913
    draw: st.DrawFn,
    *,
    root: Optional[Path] = None,
    names: st.SearchStrategy[str] | None = None,
    blacklist_entities: Optional[Container[BidsEntity | str]] = None,
    whitelist_entities: Optional[Container[BidsEntity | str]] = None,
    unique: bool = False,
    cull: bool = True,
) -> BidsDataset:
    comp1 = draw(
        bids_components(
            min_entities=2,
            max_entities=3,
            restrict_patterns=True,
            root=root,
            name=draw(names) if names is not None else None,
            blacklist_entities=blacklist_entities,
            whitelist_entities=whitelist_entities,
            unique=unique,
            cull=cull,
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


def everything_except(*excluded_types: Type[Any]) -> st.SearchStrategy[Any]:
    return (
        st.from_type(type)
        .flatmap(st.from_type)
        .filter(lambda s: not isinstance(s, excluded_types))
    )


def _is_hashable(__item: Any):
    try:
        hash(__item)
        return True
    except TypeError:
        return False


def _supports_eq(__item: Any):
    try:
        __item == 0  # type: ignore
        return True
    except Exception:
        return False


def hashables() -> st.SearchStrategy[Hashable]:
    return everything().filter(_is_hashable)


def partially_ordered() -> st.SearchStrategy[Any]:
    return everything().filter(_supports_eq)
