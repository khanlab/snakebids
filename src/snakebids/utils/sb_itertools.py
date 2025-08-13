import itertools as it
from typing import Iterable, Sequence, TypeVar

import more_itertools as itx

T = TypeVar("T")


def unpack(iterable: Iterable[T], default: Sequence[T]) -> Iterable[T]:
    """Return default if iterable has no elements.

    Allows safe unpacking of possibly empty iterables

    Parameters
    ----------
    iterable : Iterable[T]
        Iterable to unpack
    default : Sequence[T]
        Values to return if then length of iterable is 0

    Returns
    -------
    Iterable[T]
    """
    first, second = it.tee(iterable, 2)
    if itx.ilen(first):
        return second
    return default


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
