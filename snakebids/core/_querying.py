from __future__ import annotations

import abc
import functools as ft
import re
from typing import TYPE_CHECKING, Any, Final, Iterable, Mapping, Sequence, cast

import attrs
import more_itertools as itx
from bids.layout import BIDSLayout, Query
from bids.layout.models import BIDSFile
from typing_extensions import Self, TypeAlias, override

from snakebids.exceptions import ConfigError, PybidsError
from snakebids.types import FilterMap, FilterValue, InputConfig

CompiledFilter: TypeAlias = "Mapping[str, Sequence[str | Query]]"


class PostFilter:
    """Filters to apply after indexing, typically derived from the CLI.

    Currently used for supporting ``--[exclude-]participant-label``
    """

    def __init__(self):
        self.inclusions: dict[str, Sequence[str] | str] = {}
        self.exclusions: dict[str, Sequence[str] | str] = {}

    def add_filter(
        self,
        key: str,
        inclusions: Iterable[str] | str | None,
        exclusions: Iterable[str] | str | None,
    ):
        """Add entity filter based on inclusion or exclusion criteria.

        Converts a list of values to include or exclude into Pybids compatible filters.
        Exclusion filters are appropriately formatted as regex. Raises an exception if
        both include and exclude are stipulated.

        PostFilter is modified in-place.

        Parameters
        ----------
        key
            Name of entity to be filtered
        inclusions
            Values to include, values not found in this list will be excluded, by
            default ``None``
        exclusions
            Values to exclude, only values not found in this list will be included, by
            default ``None``

        Raises
        ------
        ValueError
            Raised if both include and exclude values are stipulated.
        """
        if inclusions is not None and exclusions is not None:
            msg = (
                "Cannot define both participant_label and exclude_participant_label at "
                "the same time"
            )
            raise ValueError(msg)
        if inclusions is not None:
            self.inclusions[key] = list(itx.always_iterable(inclusions))
        if exclusions is not None:
            self.exclusions[key] = self._format_exclusions(exclusions)

    def _format_exclusions(self, exclusions: Iterable[str] | str):
        # if multiple items to exclude, combine with with item1|item2|...
        exclude_string = "|".join(
            re.escape(label) for label in itx.always_iterable(exclusions)
        )
        # regex to exclude subjects
        return [f"^((?!({exclude_string})$).*)$"]


@attrs.define(slots=False)
class UnifiedFilter:
    """Manages component level and post filters."""

    component: InputConfig
    """The Component configuration defining the filters"""

    postfilters: PostFilter
    """Filters to be applied after collecting and parsing the data

    Currently only used to implement --[exclude-]participant-label, but in the future,
    may implement other such CLI args. Unlike configuration-defined filters, these
    filters apply after the dataset is indexed and queried. Thus, if a filter is set
    to an empty list, a complete, albeit empty, component may be found. This is akin to
    calling ``BidsComponent.filter`` after running ``generate_inputs``.

    For performance purposes, non-empty post-filters are applied via ``pybids.get()``
    """

    @classmethod
    def from_filter_dict(
        cls,
        filters: Mapping[str, str | bool | Sequence[str | bool]],
        postfilter: PostFilter | None = None,
    ) -> Self:
        """Patch together a UnifiedFilter based on a basic filter dict.

        Intended primarily for use in testing
        """
        wildcards: list[str] = []
        if postfilter is not None:
            wildcards.extend(postfilter.inclusions)
            wildcards.extend(postfilter.exclusions)
        return cls(
            {"filters": filters, "wildcards": wildcards}, postfilter or PostFilter()
        )

    def _has_empty_list(self, items: Iterable[Any]):
        """Check if any of the lists within iterable are empty."""
        return any(
            itx.ilen(itx.always_iterable(item, base_type=(str, dict)))  # type: ignore
            == 0
            for item in items
        )

    def _has_overlap(self, key: str):
        """Check if filter key is a wildcard and not already a prefilter."""
        return key not in self.prefilters and key in self.component.get("wildcards", [])

    @ft.cached_property
    def prefilters(self) -> FilterMap:
        """Filters defined in the component configuration and applied via pybids.

        Unlike postfilters, a prefilter set to an empty list will result in no valid
        paths found, resulting in a blank (missing) component.
        """
        filters = dict(self.component.get("filters", {}))
        # Silently remove "regex_search". This value has been blocked by a bug for the
        # since version 0.6, and even before, never fully worked properly (e.g. would
        # break if combined with --exclude-participant-label)
        if "regex_search" in filters:
            del filters["regex_search"]
        return filters

    @ft.cached_property
    def get(self) -> CompiledFilter:
        """The combination of pre- and post- filters for indexing pybids via ``.get()``.

        Includes pre-filters not annotated for regex querying and all inclusion
        post-filters. Empty post-filters are replaced with Query.ANY. This allows valid
        paths to be found and processed later. Post-filters are not applied when an
        equivalent prefilter is present

        Raises
        ------
        FilterSpecError
            When filter configuration is invalidly specified.
        """
        result = dict(_compile_filters(self.prefilters, with_regex=False))
        postfilters = self.postfilters.inclusions
        for key in self.postfilters.inclusions:
            if self._has_overlap(key):
                # if empty list filter, ensure the entity filtered is present
                result[key] = (
                    postfilters[key]
                    if itx.ilen(itx.always_iterable(postfilters[key]))
                    else [Query.ANY]
                )
        return result

    @ft.cached_property
    def search(self) -> CompiledFilter:
        """Pre-filters for indexing pybids via ``.get(regex_search=True)``.

        As with :prop:`UnifiedFilter.get`, but only prefilters labelled for regex
        matching using ``search:`` or ``match:``.

        Raises
        ------
        FilterSpecError
            When filter configuration is invalidly specified.
        """
        return dict(_compile_filters(self.prefilters, with_regex=True))

    @property
    def post_exclusions(self) -> dict[str, Sequence[str] | str]:
        """Dictionary of all post-exclusion filters."""
        return {
            key: val
            for key, val in self.postfilters.exclusions.items()
            if self._has_overlap(key)
        }

    @property
    def without_bools(self) -> Mapping[str, str | Sequence[str]]:
        """Check and typeguard to ensure filters do not contain booleans."""
        for key, val in self.get.items():
            if any(isinstance(v, Query) for v in itx.always_iterable(val)):
                msg = (
                    "Boolean filters in items with custom paths are not supported; in "
                    f"component='{key}'"
                )
                raise ValueError(msg)
        return cast("Mapping[str, str | Sequence[str]]", self.get)

    @property
    def has_empty_prefilter(self) -> bool:
        """Returns True if even one prefilter is empty."""
        return self._has_empty_list(self.prefilters.values())

    @property
    def has_empty_postfilter(self) -> bool:
        """Returns True if even one postfilter is empty."""
        return self._has_empty_list(
            filt
            for name, filt in self.postfilters.inclusions.items()
            if self._has_overlap(name)
        )


def get_matching_files(
    bids_layout: BIDSLayout,
    filters: UnifiedFilter,
) -> Iterable[BIDSFile]:
    """Query pybids layout based on provided filters.

    Supports a combination of regular and regex querying.

    Raises
    ------
    FilterSpecError
        When filter configuration is invalidly specified.
    PybidsError
        When pybids raises an error within BIDSLayout.get()
    """
    if filters.has_empty_prefilter:
        return []
    try:
        get = bids_layout.get(
            regex_search=False,
            **filters.get,
        )
        search = (
            set(bids_layout.get(regex_search=True, **filters.search))
            if filters.search
            else None
        )
    except AttributeError as err:
        msg = (
            "Pybids has encountered a problem that Snakebids cannot handle. This "
            "may indicate a missing or invalid dataset_description.json for this "
            "dataset."
        )
        raise PybidsError(msg) from err

    if search is not None:
        return [p for p in get if p in search]
    return get


@attrs.define
class FilterSpecError(Exception, abc.ABC):
    entity: str
    value: FilterValue

    requirement: Final[str] = attrs.field(
        default="Must have exactly one of {'get', 'match', search'}.", init=False
    )

    @abc.abstractmethod
    def get_config_error(self, component_name: str) -> ConfigError:
        """Return ConfigError with class-specific message."""
        ...


@attrs.define
class _TooFewKeysError(FilterSpecError):
    """Exception raised when filter specified as dictionary without keys."""

    @override
    def get_config_error(self, component_name: str) -> ConfigError:
        msg = (
            f"Filter '{self.entity}' for component '{component_name}' was specified as "
            f"a dict but was not given any keys. {self.requirement} Got: {self.value}"
        )
        return ConfigError(msg)


@attrs.define
class _TooManyKeysError(FilterSpecError):
    """Exception raised when filter specified as dictionary with more than one key."""

    @override
    def get_config_error(self, component_name: str) -> ConfigError:
        msg = (
            f"Filter '{self.entity}' for component '{component_name}' may not have "
            f"more than one key. {self.requirement} Got: {self.value}"
        )
        return ConfigError(msg)


@attrs.define
class _InvalidKeyError(FilterSpecError):
    """Exception raised when filter specified as dictionary with invalid key."""

    @override
    def get_config_error(self, component_name: str) -> ConfigError:
        msg = (
            f"Invalid query method specified for filter '{self.entity}' in component "
            f"'{component_name}'. {self.requirement} Got {self.value}"
        )
        return ConfigError(msg)


def _compile_filters(filters: FilterMap, *, with_regex: bool) -> CompiledFilter:
    """Convert filter configuration into dict consumable by pybids get().

    Raises
    ------
    FilterSpecError
        When filter configuration is invalidly specified.
    """

    def wrap(filt: str, method: str):
        if method == "get":
            return filt
        # use flag group to remove the default case-insensitivity enforced by pybids
        filt = f"(?-i:{filt})"
        if method == "match":
            # pybids only does search, so surround string filters with position anchors
            # to simulate match
            filt = f"^(?:{filt})$"
        return filt

    result: CompiledFilter = {}
    for key, raw_filter in filters.items():
        try:
            filt_type = itx.one(raw_filter.keys(), too_short=TypeError)  # type: ignore
        except ValueError as err:
            raise _TooManyKeysError(key, raw_filter) from err
        except TypeError as err:
            raise _TooFewKeysError(key, raw_filter) from err
        except AttributeError:
            if TYPE_CHECKING:
                assert not isinstance(raw_filter, dict)
            filt = raw_filter
            filt_type = "get"
        else:
            if filt_type not in {"match", "search", "get"}:
                raise _InvalidKeyError(key, raw_filter)
            if TYPE_CHECKING:
                assert isinstance(raw_filter, dict)
            filt = cast("str | bool | Sequence[str | bool]", raw_filter[filt_type])

        # these two must not be simultaneously true or false
        if with_regex == (filt_type == "get"):
            continue

        result[key] = [
            Query.ANY if f is True else Query.NONE if f is False else wrap(f, filt_type)
            for f in itx.always_iterable(filt)
        ]

    return result
