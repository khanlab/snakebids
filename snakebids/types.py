from __future__ import annotations

from enum import Enum
from pathlib import Path
from typing import Dict, Generic, Iterable, List, Mapping, Protocol, Sequence, overload

from typing_extensions import Self, TypeAlias, TypedDict, TypeVar

from snakebids.utils.containers import MultiSelectDict

_T_contra = TypeVar("_T_contra", contravariant=True)
_S_co = TypeVar("_S_co", covariant=True)


class FilterSpec(TypedDict, total=False):
    """Optional filter specification allowing regex matching."""

    get: str | bool | Sequence[str | bool]
    match: str
    search: str


FilterValue: TypeAlias = "str | bool | Sequence[str | bool] | FilterSpec"
FilterMap: TypeAlias = "Mapping[str, FilterValue]"


class InputConfig(TypedDict, total=False):
    """Configuration for a single bids component."""

    filters: FilterMap
    """Filters to pass on to :class:`BIDSLayout.get() <bids.layout.BIDSLayout>`

    Each key refers to the name of an entity. Values may take the following forms:

    * :class:`string <str>`: Restricts the entity to the exact string given
    * :class:`bool`: ``True`` requires the entity to be present (with any value).
      ``False`` requires the entity to be absent.
    * :class:`list` [:class:`str`]: List of allowable values the entity may take.

    In addition, a few special filters may be added which carry different meanings:

    * ``use_regex: True``: If present, all strings will be interpreted as regex
    * ``scope``: Restricts the scope of the component. It may take the following values:
        - ``"all"``: search everything (default behaviour)
        - ``"raw"``: only search the top-level raw dataset
        - ``"derivatives"``: only search derivative datasets
        - ``<PipelineName>``: only search derivative datasets with a matching pipeline
            name
    """

    wildcards: list[str]
    """Wildcards to allow in the component.

    Each value in the list refers to the name of an entity. If the entity is present,
    the generated :class:`~snakebids.BidsComponent` will have values of this entity
    substituted for wildcards in the :attr:`~snakebids.BidsComponent.path`, and the
    entity will be included in the :attr:`~snakebids.BidsComponent.zip_lists`.

    If the entity is not found, it will be ignored.
    """
    custom_path: str


class BinaryOperator(Protocol, Generic[_T_contra, _S_co]):
    """Callables that act on two objects of identical type."""

    def __call__(self, first: _T_contra, second: _T_contra, /) -> _S_co:  # noqa: D102
        ...


class Expandable(Protocol):
    """Protocol represents objects that hold an entity table and can expand over a path.

    Includes BidsComponent, BidsPartialComponent, and BidsComponentRow
    """

    @property
    def zip_lists(self) -> ZipList:  # noqa: D102
        ...

    def expand(  # noqa: D102
        self,
        paths: Iterable[Path | str] | Path | str,
        /,
        allow_missing: bool | str | Iterable[str] = False,
        **wildcards: str | Iterable[str],
    ) -> list[str]:
        ...

    def filter(  # noqa: D102
        self,
        *,
        regex_search: bool | str | Iterable[str] = False,
        **filters: str | Iterable[str],
    ) -> Self:
        ...


_K_contra = TypeVar("_K_contra", bound="str", contravariant=True)
_V_co = TypeVar("_V_co", covariant=True)
_Valt_co = TypeVar("_Valt_co", covariant=True)


class MultiSelectable(Protocol, Generic[_K_contra, _V_co, _Valt_co]):
    """Mappings supporting selection with multiple keys."""

    @overload
    def __getitem__(self, key: _K_contra, /) -> _V_co:
        ...

    @overload
    def __getitem__(self, key: tuple[_K_contra, ...], /) -> _Valt_co:
        ...


InputsConfig: TypeAlias = Dict[str, InputConfig]
"""Configuration for all bids components to be parsed in the app

Should be defined in the config.yaml file, by convention in a key called 'pybids_inputs'
"""

ZipList: TypeAlias = MultiSelectDict[str, List[str]]
"""Multiselectable dict mapping entity names to possible values.

All lists must be the same length. Entries in each list with the same index correspond
to the same path. Thus, the ZipList can be read like a table, where each row corresponds
to an entity, and each "column" corresponds to a path.
"""


ZipListLike: TypeAlias = Mapping[str, Sequence[str]]
"""
Generic form of a :py:data:`ZipList`

Useful for typing functions that won't mutate the ZipList or use
:class:`~snakebids.utils.utils.MultiSelectDict` capabilities. Like :class:`ZipList`,
each :class:`~typing.Sequence` must be the same length, and values in each with the same
index must correspond to the same path.
"""


class OptionalFilterType(Enum):
    """Sentinel value for CLI OPTIONAL filtering.

    This is necessary because None means no CLI filter was added.
    """

    OptionalFilter = "OptionalFilter"


OptionalFilter = OptionalFilterType.OptionalFilter
