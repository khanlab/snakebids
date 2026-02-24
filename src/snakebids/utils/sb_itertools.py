import itertools as it
from collections.abc import Iterable
from typing import TypeVar

T = TypeVar("T")


def drop(n: int, iterable: Iterable[T]) -> Iterable[T]:
    """Return all items of an iterable except the first *n* as a list.

    >>> drop(7, range(10))
    [7, 8, 9]

    If there are fewer than *n* items in the iterable, none of them are returned

    >>> drop(8, range(5))
    []
    """
    _list = list(iterable)
    return list(it.islice(_list, n, len(_list)))
