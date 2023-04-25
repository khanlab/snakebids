from __future__ import annotations

import functools as ft
import importlib.resources
import json
import operator as op
import re
from pathlib import Path
from typing import Any, Callable, Iterable, Mapping, Sequence, TypeVar, cast, overload

import attrs
import more_itertools as itx
from typing_extensions import Protocol, Self

from snakebids import types
from snakebids.utils.user_property import UserProperty

T = TypeVar("T")


@ft.lru_cache(None)
def read_bids_tags(bids_json: Path | None = None) -> dict[str, dict[str, str]]:
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
            bids_tags = json.load(infile)
        return bids_tags
    with importlib.resources.open_text("snakebids.resources", "bids_tags.json") as file:
        bids_tags = json.load(file)
    return bids_tags


@attrs.frozen(hash=True)
class BidsEntity:
    """Bids entities with tag and wildcard representations"""

    entity: str = attrs.field(converter=str)

    def __str__(self) -> str:
        return self.entity

    def __eq__(self, other: Any) -> bool:
        if isinstance(other, str):
            return self.entity == other
        if isinstance(other, BidsEntity):
            return self.entity == other.entity
        return False

    @property
    def tag(self) -> str:
        """Get the bids tag version of the entity

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
        """Get regex of acceptable value matches

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
        """regex str to search before value in paths"""
        tags = read_bids_tags()
        return (
            tags[self.entity]["before"]
            if self.entity in tags and "before" in tags[self.entity]
            else f"{self.tag}-"
        )

    @property
    def after(self) -> str:
        """regex str to search after value in paths"""
        tags = read_bids_tags()
        return (
            tags[self.entity]["after"]
            if self.entity in tags and "after" in tags[self.entity]
            else ""
        )

    @property
    def regex(self) -> re.Pattern[str]:
        """Complete pattern to match when searching in paths

        Contains three capture groups, the first corresponding to "before", the second
        to "value", and the third to "after"
        """
        return re.compile(f"({self.before})({self.match})({self.after})")

    @property
    def wildcard(self) -> str:
        """Get the snakebids {wildcard}

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
        """Return the entity associated with the given tag, if found

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


def matches_any(
    item: T,
    match_list: Iterable[T],
    match_func: types.BinaryOperator[T, object],
    *args: Any,
) -> bool:
    for match in match_list:
        if match_func(match, item, *args):
            return True
    return False


def get_match_search_func(
    match_list: Iterable[Any], match_func: Callable[[Any, Any], Any]
) -> Callable[[Any], bool]:
    """Return a match search function suitable for use in filter

    Parameters
    ----------
    match_list : list
        list of items to search for matches
    match_func : callable
        Any callable that takes two args and returns truthy/falsy values. The first arg
        will be drawn from match_list, the second will be the value being matched

    Returns
    -------
    callable
        Takes as a single arg a value to be matched. It will be compared to every item
        in match list using match_func
    """
    match_list = list(match_list)

    def inner(item: Any):
        return matches_any(item, match_list, match_func)

    return inner


class BidsParseError(Exception):
    """Exception raised for errors encountered in the parsing of Bids paths"""

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
) -> Callable[[Callable[[Any], T]], "UserProperty[T]"]:
    """Set property as an alias for another property

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

    def inner(__func: Callable[[Any], T]) -> "UserProperty[T]":
        alias = UserProperty(__func)
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


def surround(__s: Iterable[str] | str, __object: str) -> Iterable[str]:
    """Surround a string or each string in an iterable with characters"""
    for item in itx.always_iterable(__s):
        yield __object + item + __object


_K = TypeVar("_K", bound="str")
_V = TypeVar("_V")


class MultiSelectDict(types.UserDictPy37[_K, _V]):
    """Dict supporting selection of multiple keys using tuples

    If a single key is given, the item associated with that key is returned just as in a
    regular dict. If multiple, comma-seperated keys are given, (e.g. a tuple), a new
    ``MultiSelectDict`` will be returned containing the keys given and their values:

        >>> mydict = MultiSelectDict({
        ...     "foo": "bar",
        ...     "hello": "world",
        ...     "fee": "fie",
        ... })
        >>> mydict["foo"]
        'bar'
        >>> mydict["foo", "hello"]
        {'foo': 'bar', 'hello': 'world'}

    The new ``MultiSelectDict`` is a "view" of the original data. Any mutations made to
    the values will be reflected in the original object:

        >>> mydict = MultiSelectDict({
        ...     "foo": [1, 2, 3, 4],
        ...     "hello": "world",
        ...     "fee": "fie",
        ... })
        >>> view = mydict["foo", "hello"]
        >>> view["foo"].append(5)
        >>> mydict
        {'foo': [1, 2, 3, 4, 5], 'hello': 'world', 'fee': 'fie'}

    The keys of the new ``MultiSelectDict`` will be inserted in the order provided in
    the selector

        >>> mydict = MultiSelectDict({
        ...     "foo": "bar",
        ...     "hello": "world",
        ...     "fee": "fie",
        ... })
        >>> mydict["hello", "foo"]
        {'hello': 'world', 'foo': 'bar'}

    Tuples of length 1 will still return a ``MultiSelectDict``:

        >>> mydict = MultiSelectDict({
        ...     "foo": [1, 2, 3, 4],
        ...     "hello": "world",
        ...     "fee": "fie",
        ... })
        >>> mydict["foo",]
        {'foo': [1, 2, 3, 4]}

    If a key is given multiple times, the extra instances of the key will be ignored:

        >>> mydict = MultiSelectDict({
        ...     "foo": [1, 2, 3, 4],
        ...     "hello": "world",
        ...     "fee": "fie",
        ... })
        >>> mydict["foo", "foo", "foo"]
        {'foo': [1, 2, 3, 4]}
    """

    @overload
    def __getitem__(self, __key: _K) -> _V:
        ...

    @overload
    def __getitem__(self, __key: tuple[_K]) -> Self:
        ...

    def __getitem__(self, __key: _K | tuple[_K]) -> _V | Self:
        if isinstance(__key, tuple):
            # Use dict.fromkeys for de-duplication
            return self.__class__({key: self[key] for key in dict.fromkeys(__key)})
        return super().__getitem__(__key)


def zip_list_eq(__first: types.ZipListLike, __second: types.ZipListLike):
    """Compare two zip lists, allowing the order of columns to be irrelevant"""

    def sorted_items(dictionary: Mapping[str, Sequence[str]]):
        return sorted(dictionary.items(), key=op.itemgetter(0))

    def get_values(zlist: types.ZipListLike):
        return cast("tuple[list[str]]", list(zip(*sorted_items(zlist)))[1])

    if not __first and not __second:
        return True

    if set(__first) != set(__second):
        return False

    first_items = get_values(__first)
    second_items = get_values(__second)

    return set(zip(*first_items)) == set(zip(*second_items))
