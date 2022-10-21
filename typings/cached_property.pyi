from typing import TypeVar

from snakebids.utils.sb_typing import UserProperty

T = TypeVar("T")

class cached_property(UserProperty[T]): ...
