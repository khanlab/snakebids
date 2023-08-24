from __future__ import annotations

import functools as ft
import json
import operator as op
import os
import re
import sys
from os import PathLike
from pathlib import Path
from typing import (
    Any,
    Callable,
    Generic,
    Iterable,
    Iterator,
    Mapping,
    Sequence,
    TypeVar,
    cast,
    overload,
)

import attrs
import importlib_resources as impr
import more_itertools as itx
from typing_extensions import (
    NotRequired,
    Protocol,
    Self,
    SupportsIndex,
    TypeAlias,
    TypedDict,
)

from snakebids import resources, types
from snakebids.utils.user_property import UserProperty

_T = TypeVar("_T")


class BidsTag(TypedDict):
    tag: str
    before: NotRequired[str]
    match: str
    after: NotRequired[str]
    leader: NotRequired[bool]


BidsTags: TypeAlias = "dict[str, BidsTag]"


DEPRECATION_FLAG = "<!DEPRECATED!>"
"""Sentinel string to mark deprecated config features"""


@ft.lru_cache(None)
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
            bids_tags = json.load(infile)
        return bids_tags
    return json.loads(impr.files(resources).joinpath("bids_tags.json").read_text())


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

    def __lt__(self, other: BidsEntity | str):
        if isinstance(other, str):
            return self.entity < other
        return self.entity < other.entity

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
        return tags.get(self.entity, {}).get("before", f"{self.tag}-")

    @property
    def after(self) -> str:
        """regex str to search after value in paths"""
        tags = read_bids_tags()
        return tags.get(self.entity, {}).get("after", "")

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

    @classmethod
    def normalize(cls, __item: str | BidsEntity) -> BidsEntity:
        """Return the entity associated with the given item, if found

        Supports both strings and BidsEntities as input. Unlike the constructor, if a
        tag name is given, the associated entity will be returned. If no associated
        entity is found, the tag itself is used as the entity name

        Parameters
        ----------
        tag : str
            tag to search

        Returns
        -------
        BidsEntity
        """
        if isinstance(__item, BidsEntity):
            return __item
        return cls.from_tag(__item)


def matches_any(
    item: _T,
    match_list: Iterable[_T],
    match_func: types.BinaryOperator[_T, object],
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
) -> Callable[[Callable[[Any], _T]], "UserProperty[_T]"]:
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

    def inner(__func: Callable[[Any], _T]) -> "UserProperty[_T]":
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

    return sorted(zip(*first_items)) == sorted(zip(*second_items))


def get_first_dir(path: str) -> str:
    """Return the top level directory in a path

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
    return Path(path).resolve()


_T_co = TypeVar("_T_co", covariant=True)


class ImmutableList(Sequence[_T_co], Generic[_T_co]):
    """Subclassable tuple equivalent

    Mimics a tuple in every way, but readily supports subclassing. Data is stored on a
    private attribute ``_data``. Subclasses must not override this attribute. To avoid
    accidental modification, subclasses should avoid interacting with ``_data``, using
    the relevant ``super()`` calls to access internal data instead (e.g. use
    ``super().__getitem__(index)`` rather than ``self._data[index]``).

    Unlike tuples, only a single type parameter is supported. In other words,
    ``ImmutableList`` cannot be specified via type hints as a fixed length sequence
    containing heterogenous items. A tuple specified as ``tuple[str, int, str]`` would
    be specified as ``ImmutableList[str | int]``
    """

    def __init__(self, __iterable: Iterable[_T_co] = tuple()):
        self._data = tuple(__iterable)

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}({list(self._data)})"

    def __contains__(self, __item: object) -> bool:
        return __item in self._data

    def __hash__(self):
        return hash(self._data)

    def __iter__(self) -> Iterator[_T_co]:
        return iter(self._data)

    def __reversed__(self) -> Iterator[_T_co]:
        return reversed(self._data)

    def __len__(self) -> int:
        return len(self._data)

    def __bool__(self) -> bool:
        return bool(self._data)

    @overload
    def __getitem__(self, __key: SupportsIndex) -> _T_co:
        ...

    @overload
    def __getitem__(self, __key: slice) -> Self:
        ...

    def __getitem__(self, __item: SupportsIndex | slice) -> _T_co | Self:
        if isinstance(__item, slice):
            return self.__class__(self._data[__item])
        return self._data[__item]

    def __lt__(self, __value: tuple[_T_co, ...] | Self) -> bool:
        if isinstance(__value, tuple):
            return self._data < __value
        if isinstance(__value, ImmutableList):
            return self._data < __value._data
        return False

    def __le__(self, __value: tuple[_T_co, ...] | Self) -> bool:
        if isinstance(__value, tuple):
            return self._data <= __value
        if isinstance(__value, ImmutableList):
            return self._data <= __value._data
        return False

    def __gt__(self, __value: tuple[_T_co, ...] | Self) -> bool:
        if isinstance(__value, tuple):
            return self._data > __value
        if isinstance(__value, ImmutableList):
            return self._data > __value._data
        return False

    def __ge__(self, __value: tuple[_T_co, ...] | Self) -> bool:
        if isinstance(__value, tuple):
            return self._data >= __value
        if isinstance(__value, ImmutableList):
            return self._data >= __value._data
        return False

    def __eq__(self, __value: object) -> bool:
        if isinstance(__value, tuple):
            return self._data == __value
        if isinstance(__value, ImmutableList):
            return self._data == __value._data  # type: ignore
        return False

    def __add__(self, __value: tuple[_T_co, ...] | ImmutableList[_T_co]) -> Self:
        if isinstance(__value, ImmutableList):
            return self.__class__(self._data + __value._data)
        return self.__class__(self._data + __value)

    def __mul__(self, __value: SupportsIndex) -> Self:
        return self.__class__(self._data * __value)

    def __rmul__(self, __value: SupportsIndex) -> Self:
        return self.__class__(__value * self._data)

    def count(self, value: Any) -> int:
        return self._data.count(value)

    def index(
        self, value: Any, start: SupportsIndex = 0, stop: SupportsIndex = sys.maxsize
    ) -> int:
        return self._data.index(value, start, stop)


def get_wildcard_dict(__entities: str | Iterable[str]) -> dict[str, str]:
    """Turn entity strings into wildcard dicts as {"entity": "{entity}"}"""
    return {entity: f"{{{entity}}}" for entity in itx.always_iterable(__entities)}
