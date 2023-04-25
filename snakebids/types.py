from __future__ import annotations

from collections.abc import Hashable
from typing import Dict, Generic, List, Mapping, Sequence

from typing_extensions import TYPE_CHECKING, Protocol, TypeAlias, TypedDict, TypeVar

_T_contra = TypeVar("_T_contra", contravariant=True)
_S_co = TypeVar("_S_co", covariant=True)


class InputConfig(TypedDict, total=False):
    """Configuration for a single bids component"""

    filters: dict[str, str | bool | list[str]]
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
    def __call__(self, __first: _T_contra, __second: _T_contra) -> _S_co:
        ...


# Hack to make userdicts subscriptable in python 3.7. Can remove when we drop support
# for that version
_K = TypeVar("_K", bound=Hashable)
_V = TypeVar("_V")
if TYPE_CHECKING:

    class UserDictPy37(Dict[_K, _V]):
        pass

else:

    class UserDictPy37(dict, Generic[_K, _V]):
        pass


# for py37, we need to import this AFTER we initialize UserDictPy37 to avoid a circular
# import
from snakebids.utils import utils  # noqa: E402

InputsConfig: TypeAlias = Dict[str, InputConfig]
"""Configuration for all bids components to be parsed in the app

Should be defined in the config.yaml file, by convention in a key called 'pybids_inputs'
"""

ZipList: TypeAlias = "utils.MultiSelectDict[str, List[str]]"
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
