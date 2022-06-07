from typing import Callable, TypeVar

T = TypeVar("T")

def cached_property(func: Callable[..., T]) -> T: ...
