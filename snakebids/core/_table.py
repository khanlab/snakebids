from __future__ import annotations

from typing import TYPE_CHECKING, Iterable, Mapping

import attrs
import more_itertools as itx

from snakebids.io.printing import format_zip_lists
from snakebids.types import ZipList, ZipListLike
from snakebids.utils.containers import ContainerBag, MultiSelectDict, RegexContainer

if TYPE_CHECKING:

    def wcard_tuple(x: Iterable[str]) -> tuple[str, ...]:
        return tuple(x)

    def entries_list(x: Iterable[tuple[str, ...]]) -> list[tuple[str, ...]]:
        return list(x)

    def liststr() -> list[str]:
        return []
else:
    wcard_tuple = tuple
    entries_list = list
    liststr = list


@attrs.frozen(kw_only=True)
class BidsTable:
    """Container holding the entries of a BidsComponent."""

    wildcards: tuple[str, ...] = attrs.field(converter=wcard_tuple)
    entries: list[tuple[str, ...]] = attrs.field(converter=entries_list)

    def __bool__(self):
        """Return True if one or more entries, otherwise False."""
        return bool(self.entries)

    def __eq__(self, other: BidsTable | object):
        if not isinstance(other, self.__class__):
            return False
        if set(self.wildcards) != set(other.wildcards):
            return False
        if len(self.entries) != len(other.entries):
            return False
        if self.wildcards == other.wildcards:
            return sorted(self.entries) == sorted(other.entries)
        ixs = [other.wildcards.index(w) for w in self.wildcards]
        entries = self.entries.copy()
        try:
            for entry in other.entries:
                sorted_entry = tuple(entry[i] for i in ixs)
                entries.remove(sorted_entry)
        except ValueError:
            return False
        return True

    @classmethod
    def from_dict(cls, d: ZipListLike):
        """Construct BidsTable from a mapping of entities to value lists."""
        lengths = {len(val) for val in d.values()}
        if len(lengths) > 1:
            msg = "each entity must have the same number of values"
            raise ValueError(msg)
        return cls(wildcards=d.keys(), entries=zip(*d.values()))

    def to_dict(self) -> ZipList:
        """Convert into a zip_list."""
        if not self.entries:
            return MultiSelectDict(zip(self.wildcards, itx.repeatfunc(liststr)))
        return MultiSelectDict(zip(self.wildcards, map(list, zip(*self.entries))))

    def pformat(self, max_width: int | float | None = None, tabstop: int = 4) -> str:
        """Pretty-format."""
        return format_zip_lists(self.to_dict(), max_width=max_width, tabstop=tabstop)

    def get(self, wildcard: str):
        """Get values for a single wildcard."""
        index = self.wildcards.index(wildcard)
        return [entry[index] for entry in self.entries]

    def pick(self, wildcards: Iterable[str]):
        """Select wildcards without deduplication."""
        # Use dict.fromkeys for de-duplication to preserve order
        unique_keys = list(dict.fromkeys(wildcards))
        indices = [self.wildcards.index(w) for w in unique_keys]

        entries = [tuple(entry[i] for i in indices) for entry in self.entries]

        return self.__class__(wildcards=unique_keys, entries=entries)

    def filter(
        self,
        filters: Mapping[str, Iterable[str] | str],
        regex_search: bool = False,
    ):
        """Apply filtering operation."""
        valid_filters = set(self.wildcards)
        if regex_search:
            filter_sets = {
                self.wildcards.index(key): ContainerBag(
                    *(RegexContainer(r) for r in itx.always_iterable(vals))
                )
                for key, vals in filters.items()
                if key in valid_filters
            }
        else:
            filter_sets = {
                self.wildcards.index(key): set(itx.always_iterable(vals))
                for key, vals in filters.items()
                if key in valid_filters
            }

        keep = [
            entry
            for entry in self.entries
            if all(
                i not in filter_sets or val in filter_sets[i]
                for i, val in enumerate(entry)
            )
        ]

        return self.__class__(wildcards=self.wildcards, entries=keep)
