import itertools as it
from typing import Iterable, Sequence, TypeVar

import more_itertools as itx

T = TypeVar("T")


def unpack(iterable: Iterable[T], default: Sequence[T]) -> Iterable[T]:
    """Return default if iterable has no elements

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


# pylint: disable=invalid-name
def drop(n: int, iterable: Iterable[T]) -> Iterable[T]:
    first, second = it.tee(iterable, 2)
    keep = itx.ilen(first) - n
    return itx.take(max(0, keep), second)
