from __future__ import annotations

import functools as ft
import json
import operator as op
import os
import re
from os import PathLike
from pathlib import Path
from typing import Any, Callable, Iterable, Mapping, Protocol, Sequence, TypeVar, cast

import attrs
import importlib_resources as impr
import more_itertools as itx
from typing_extensions import NotRequired, TypeAlias, TypedDict

from snakebids import resources, types
from snakebids.utils.user_property import UserProperty

_T = TypeVar("_T")


class BidsTag(TypedDict):
    """Interface for BidsTag configuration."""

    tag: str
    before: NotRequired[str]
    match: str
    after: NotRequired[str]
    leader: NotRequired[bool]


BidsTags: TypeAlias = "dict[str, BidsTag]"


DEPRECATION_FLAG = "<!DEPRECATED!>"
"""Sentinel string to mark deprecated config features"""


@ft.lru_cache
def read_bids_tags(bids_json: Path | None = None) -> BidsTags:
    """Read the bids tags we are aware of from a JSON file.

    This is used specifically for compatibility with pybids, since some tag keys
    are different from how they appear in the file name, e.g. ``subject`` for
    ``sub``, and ``acquisition`` for ``acq``.

    Parameters
    ----------
    bids_json
        Path to JSON file to use, if not specified will use
        ``bids_tags.json`` in the snakebids module.

    Returns
    -------
    dict:
        Dictionary of bids tags
    """
    if bids_json:
        with bids_json.open("r") as infile:
            return json.load(infile)
    return json.loads(impr.files(resources).joinpath("bids_tags.json").read_text())


@attrs.frozen(hash=True)
class BidsEntity:
    """Bids entities with tag and wildcard representations."""

    entity: str = attrs.field(converter=str)

    def __str__(self) -> str:
        return self.entity

    def __eq__(self, other: Any) -> bool:
        if isinstance(other, str):
            return self.entity == other
        if isinstance(other, BidsEntity):
            return self.entity == other.entity
        return False

    def __lt__(self, other: BidsEntity | str):
        if isinstance(other, str):
            return self.entity < other
        return self.entity < other.entity

    @property
    def tag(self) -> str:
        """Get the bids tag version of the entity.

        For entities in the bids spec, the tag is the short version of the entity
        name. Otherwise, the tag is equal to the entity.
        """
        tags = read_bids_tags()
        return (
            tags[self.entity]["tag"]
            if self.entity in tags and "tag" in tags[self.entity]
            else self.entity
        )

    @property
    def match(self) -> str:
        """Get regex of acceptable value matches.

        If no pattern is associated with the entity, the default pattern is a word with
        letters and numbers
        """
        tags = read_bids_tags()
        return (
            tags[self.entity]["match"]
            if self.entity in tags and "match" in tags[self.entity]
            else "[a-zA-Z0-9]+"
        )

    @property
    def before(self) -> str:
        """Regex str to search before value in paths."""
        tags = read_bids_tags()
        # Need to explicitly annotate the default here and in .after because tags is a
        # dict of `BidsTag`, a `TypedDict`. So putting unannotated dicts directly as
        # default leads to a type error
        _def: dict[Any, Any] = {}
        return tags.get(self.entity, _def).get("before", f"{self.tag}-")

    @property
    def after(self) -> str:
        """Regex str to search after value in paths."""
        tags = read_bids_tags()
        # See note in .before
        _def: dict[Any, Any] = {}
        return tags.get(self.entity, _def).get("after", "")

    @property
    def regex(self) -> re.Pattern[str]:
        """Complete pattern to match when searching in paths.

        Contains three capture groups, the first corresponding to "before", the second
        to "value", and the third to "after"
        """
        return re.compile(f"({self.before})({self.match})({self.after})")

    @property
    def wildcard(self) -> str:
        """Get the snakebids {wildcard}.

        The wildcard is generally equal to the tag, i.e. the short version of the entity
        name, except for subject and session, which use the full name name. This is to
        ensure compatibility with the bids function
        """
        # HACK FIX FOR acq vs acquisition etc -- should
        # eventually update the bids() function to also use
        # bids_tags.json, where e.g. acquisition -> acq is
        # defined.. -- then, can use wildcard_name instead
        # of out_name..
        if self.entity in ["subject", "session"]:
            return self.entity
        return self.tag

    @classmethod
    def from_tag(cls, tag: str) -> BidsEntity:
        """Return the entity associated with the given tag, if found.

        If not associated entity is found, the tag itself is used as the entity name

        Parameters
        ----------
        tag : str
            tag to search

        Returns
        -------
        BidsEntity
        """
        for entity, props in read_bids_tags().items():
            if props.get("tag", None) == tag:
                return cls(entity)
        return cls(tag)

    @classmethod
    def normalize(cls, item: str | BidsEntity, /) -> BidsEntity:
        """Return the entity associated with the given item, if found.

        Supports both strings and BidsEntities as input. Unlike the constructor, if a
        tag name is given, the associated entity will be returned. If no associated
        entity is found, the tag itself is used as the entity name

        Parameters
        ----------
        item
            tag to search
        """
        if isinstance(item, BidsEntity):
            return item
        return cls.from_tag(item)


def matches_any(
    item: _T,
    match_list: Iterable[_T],
    match_func: types.BinaryOperator[_T, object],
    *args: Any,
) -> bool:
    """Test if item matches any of the items in match_list.

    Parameters
    ----------
    item
        Item to test
    match_list
        Items to compare with
    match_func
        Function to test equality. Defaults to basic equality (``==``) check
    """
    return any(match_func(match, item, *args) for match in match_list)


class BidsParseError(Exception):
    """Exception raised for errors encountered in the parsing of Bids paths."""

    def __init__(self, path: str, entity: BidsEntity) -> None:
        self.path = path
        self.entity = entity
        super().__init__(path, entity)


class _Documented(Protocol):
    __doc__: str


def property_alias(
    prop: _Documented,
    label: str | None = None,
    ref: str | None = None,
    copy_extended_docstring: bool = False,
) -> Callable[[Callable[[Any], _T]], UserProperty[_T]]:
    """Set property as an alias for another property.

    Copies the docstring from the aliased property to the alias

    Parameters
    ----------
    prop : property
        Property to alias
    label
        Text to use in link to aliased property
    ref
        Name of the property to alias
    copy_extended_docstring
        If True, copies over the entire docstring, in addition to the summary line

    Returns
    -------
    property
    """

    def inner(func: Callable[[Any], _T], /) -> UserProperty[_T]:
        alias = UserProperty(func)
        if label:
            link = f":attr:`{label} <{ref}>`" if ref else label
        else:
            link = f":attr:`{ref}`" if ref else None
        labeltxt = f"Alias of {link}\n\n" if link else ""
        if copy_extended_docstring:
            alias.__doc__ = f"{labeltxt}{prop.__doc__}"
        else:
            alias.__doc__ = f"{labeltxt}{prop.__doc__.splitlines()[0]}"
        return alias

    return inner


def surround(s: Iterable[str] | str, object_: str, /) -> Iterable[str]:
    """Surround a string or each string in an iterable with characters."""
    for item in itx.always_iterable(s):
        yield object_ + item + object_


def zip_list_eq(first: types.ZipListLike, second: types.ZipListLike, /):
    """Compare two zip lists, allowing the order of columns to be irrelevant."""

    def sorted_items(dictionary: Mapping[str, Sequence[str]]):
        return sorted(dictionary.items(), key=op.itemgetter(0))

    def get_values(zlist: types.ZipListLike):
        return cast("tuple[list[str]]", list(zip(*sorted_items(zlist)))[1])

    if not first and not second:
        return True

    if set(first) != set(second):
        return False

    first_items = get_values(first)
    second_items = get_values(second)

    return sorted(zip(*first_items)) == sorted(zip(*second_items))


def get_first_dir(path: str) -> str:
    """Return the top level directory in a path.

    If absolute, return the root. This function is necessary to handle paths with
    ``./``, as ``pathlib.Path`` filters this out.
    """
    if os.path.isabs(path):
        return Path(path).root
    parent, child = os.path.split(path)
    if parent:
        return get_first_dir(parent)
    return child


def to_resolved_path(path: str | PathLike[str]):
    """Convert provided object into resolved path."""
    return Path(path).resolve()


def get_wildcard_dict(entities: str | Iterable[str], /) -> dict[str, str]:
    """Turn entity strings into wildcard dicts as {"entity": "{entity}"}."""
    return {entity: f"{{{entity}}}" for entity in itx.always_iterable(entities)}


def entity_to_wildcard(entities: str | Iterable[str], /):
    """Turn entity strings into wildcard dicts as {"entity": "{entity}"}."""
    return {entity: f"{{{entity}}}" for entity in itx.always_iterable(entities)}
