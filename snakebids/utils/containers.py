from __future__ import annotations

import re
import sys
from typing import (
    TYPE_CHECKING,
    Any,
    AnyStr,
    Container,
    Dict,
    Generic,
    Hashable,
    Iterable,
    Iterator,
    Sequence,
    SupportsIndex,
    TypeVar,
    overload,
)

from typing_extensions import Self, override

_T = TypeVar("_T")
_K = TypeVar("_K", bound="str")
_T_co = TypeVar("_T_co", covariant=True)

# Hack to make dicts subscriptable in python 3.8. Can remove when we drop support
# for that version
_H = TypeVar("_H", bound=Hashable)
_V = TypeVar("_V")
if TYPE_CHECKING:

    class UserDictPy38(Dict[_H, _V]):
        """Wrapper around dict, used for subclassing with static typing."""

else:

    class UserDictPy38(dict, Generic[_K, _V]):
        """Wrapper around dict, used for subclassing with static typing."""


class ImmutableList(Sequence[_T_co], Generic[_T_co]):
    """Subclassable tuple equivalent.

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

    def __init__(self, iterable: Iterable[_T_co] = (), /):
        self._data = tuple(iterable)

    @override
    def __repr__(self) -> str:
        return f"{self.__class__.__name__}({list(self._data)})"

    @override
    def __contains__(self, item: object, /) -> bool:
        return item in self._data

    @override
    def __hash__(self):
        return hash(self._data)

    @override
    def __iter__(self) -> Iterator[_T_co]:
        return iter(self._data)

    @override
    def __reversed__(self) -> Iterator[_T_co]:
        return reversed(self._data)

    @override
    def __len__(self) -> int:
        return len(self._data)

    def __bool__(self) -> bool:
        return bool(self._data)

    @overload
    def __getitem__(self, index: int) -> _T_co:
        ...

    @overload
    def __getitem__(self, index: slice) -> Self:
        ...

    @override
    def __getitem__(self, index: int | slice) -> _T_co | Self:
        if isinstance(index, slice):
            return self.__class__(self._data[index])
        return self._data[index]

    def __lt__(self, value: tuple[_T_co, ...] | Self, /) -> bool:
        if isinstance(value, ImmutableList):
            return self._data < value._data
        return self._data < value

    def __le__(self, value: tuple[_T_co, ...] | Self, /) -> bool:
        if isinstance(value, ImmutableList):
            return self._data <= value._data
        return self._data <= value

    def __gt__(self, value: tuple[_T_co, ...] | Self, /) -> bool:
        if isinstance(value, ImmutableList):
            return self._data > value._data
        return self._data > value

    def __ge__(self, value: tuple[_T_co, ...] | Self, /) -> bool:
        if isinstance(value, ImmutableList):
            return self._data >= value._data
        return self._data >= value

    @override
    def __eq__(self, value: object, /) -> bool:
        if isinstance(value, tuple):
            return self._data == value
        if isinstance(value, ImmutableList):
            return self._data == value._data  # type: ignore
        return False

    def __add__(self, value: tuple[_T_co, ...] | ImmutableList[_T_co], /) -> Self:
        if isinstance(value, ImmutableList):
            return self.__class__(self._data + value._data)
        return self.__class__(self._data + value)

    def __mul__(self, value: SupportsIndex, /) -> Self:
        return self.__class__(self._data * value)

    def __rmul__(self, value: SupportsIndex, /) -> Self:
        return self.__class__(value * self._data)

    @override
    def count(self, value: Any) -> int:
        return self._data.count(value)

    @override
    def index(
        self, value: Any, start: SupportsIndex = 0, stop: SupportsIndex = sys.maxsize
    ) -> int:
        return self._data.index(value, start, stop)


class RegexContainer(Generic[AnyStr], Container[AnyStr]):
    """Container that tests if a string matches a regex using the ``in`` operator.

    Constructed with a regex expression. Supports inclusion tests for strings using
    ``in``. Strings matching the regex (using ``re.match``) will return ``True``
    """

    def __init__(self, template: AnyStr | re.Pattern[AnyStr]):
        self.template: re.Pattern[AnyStr] = re.compile(template)

    def __contains__(self, x: object, /):
        if isinstance(x, type(self.template.pattern)):
            return self.template.match(x) is not None
        return False


class ContainerBag(Container[_T]):
    """Container to hold other containers.

    Useful because list(Container) isn't guaranteed to work, so this lets us merge
    Containers in a type safe way.
    """

    def __init__(self, *entries: Container[_T]):
        self.entries = entries

    def __contains__(self, x: object, /) -> bool:
        return any(x in entry for entry in self.entries)


class MultiSelectDict(UserDictPy38[_K, _T]):
    """Dict supporting selection of multiple keys using tuples.

    If a single key is given, the item associated with that key is returned just as in a
    regular dict. If multiple, comma-seperated keys are given, (e.g. a tuple), a new
    ``MultiSelectDict`` will be returned containing the keys given and their values::

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
    the values will be reflected in the original object::

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
    the selector::

        >>> mydict = MultiSelectDict({
        ...     "foo": "bar",
        ...     "hello": "world",
        ...     "fee": "fie",
        ... })
        >>> mydict["hello", "foo"]
        {'hello': 'world', 'foo': 'bar'}

    Tuples of length 1 will still return a ``MultiSelectDict``::

        >>> mydict = MultiSelectDict({
        ...     "foo": [1, 2, 3, 4],
        ...     "hello": "world",
        ...     "fee": "fie",
        ... })
        >>> mydict["foo",]
        {'foo': [1, 2, 3, 4]}

    If a key is given multiple times, the extra instances of the key will be ignored::

        >>> mydict = MultiSelectDict({
        ...     "foo": [1, 2, 3, 4],
        ...     "hello": "world",
        ...     "fee": "fie",
        ... })
        >>> mydict["foo", "foo", "foo"]
        {'foo': [1, 2, 3, 4]}
    """

    @overload
    def __getitem__(self, key: _K, /) -> _T:
        ...

    @overload
    def __getitem__(self, key: tuple[_K, ...], /) -> Self:
        ...

    def __getitem__(self, key: _K | tuple[_K, ...], /) -> _T | Self:
        if isinstance(key, tuple):
            # Use dict.fromkeys for de-duplication
            return self.__class__({key: self[key] for key in dict.fromkeys(key)})
        return super().__getitem__(key)
