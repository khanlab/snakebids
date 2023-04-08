# pylint: disable-all
from __future__ import annotations

from collections.abc import Hashable
from typing import Generic

from typing_extensions import TYPE_CHECKING, TypeAlias, TypedDict, TypeVar

if TYPE_CHECKING:
    # This TYPE_CHECKING guard can be removed when py37 support is dropped, as we won't
    # have a circular import anymore
    from snakebids.utils import utils


class InputConfig(TypedDict, total=False):
    """Configuration passed in snakebids.yaml file"""

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


InputsConfig: TypeAlias = "dict[str, InputConfig]"

ZipLists: TypeAlias = "utils.MultiSelectDict[str, list[str]]"

# Hack to make userdicts subscriptable in python 3.7. Can remove when we drop support
# for that version
_K = TypeVar("_K", bound=Hashable)
_V = TypeVar("_V")
if TYPE_CHECKING:

    class UserDictPy37(dict[_K, _V]):
        pass

else:

    class UserDictPy37(dict, Generic[_K, _V]):
        pass
