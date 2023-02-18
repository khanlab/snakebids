from datetime import date
from typing import Callable, Generic, ParamSpec, TypeVar

P = ParamSpec("P")
T = TypeVar("T", bound=Callable)

class DeprecatedWarning(DeprecationWarning, Generic[T]):
    def __init__(
        self,
        function: T,
        deprecated_in: str,
        removed_in: str | date | None,
        details: str = ...,
    ): ...
    def __str__(self) -> str: ...

class UnsupportedWarning(DeprecatedWarning):
    def __str__(self) -> str: ...

def deprecated(
    deprecated_in: str | None = ...,
    removed_in: str | date | None = ...,
    current_version: str | None = ...,
    details: str | None = ...,
    admonition: str | None = ...,
) -> Callable[[T], T]: ...
def fail_if_not_removed(method: T) -> T: ...
