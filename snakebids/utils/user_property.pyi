from typing import Any, Callable, Generic, TypeVar

T = TypeVar("T")

class UserProperty(Generic[T]):
    def __init__(self, method: Callable[[Any], T], /) -> None: ...
    def __get__(self, obj: Any, objtype: type[Any] = ...) -> T: ...
