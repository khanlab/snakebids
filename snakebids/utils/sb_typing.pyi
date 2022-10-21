from typing import Any, Callable, Generic, Type, TypeVar

T = TypeVar("T")

class UserProperty(Generic[T]):
    def __init__(self, __method: Callable[[Any], T]): ...
    def __get__(self, obj: Any, objtype: Type[Any] = ...) -> T: ...
