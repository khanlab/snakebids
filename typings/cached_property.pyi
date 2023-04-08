from typing import TypeVar

from snakebids.utils.user_property import UserProperty

T = TypeVar("T")

class cached_property(UserProperty[T]): ...
