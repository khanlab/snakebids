from __future__ import annotations

import argparse
import sys
from typing import (
    TYPE_CHECKING,
    Any,
    Callable,
    Generic,
    Iterable,
    TypedDict,
    TypeVar,
    overload,
)

from typing_extensions import Required, TypeAlias, Unpack

if TYPE_CHECKING:
    from argparse import _SUPPRESS_T, _ActionStr, _ArgumentGroup, _NArgsStr

_T = TypeVar("_T")


if sys.version_info >= (3, 11):

    class AddArgumentArgs(TypedDict, Generic[_T], total=False):
        """Arguments for add_argument, with ``dest`` required."""

        action: _ActionStr | type[argparse.Action]
        nargs: int | _NArgsStr | _SUPPRESS_T | None
        const: Any
        default: Any
        type: Callable[[Any], _T] | argparse.FileType
        choices: Iterable[_T] | None
        required: bool
        help: str | None
        metavar: str | tuple[str, ...] | None
        dest: Required[str]
        version: str

    AnyArgumentArgs: TypeAlias = "AddArgumentArgs[Any]"

else:

    class AddArgumentArgs(TypedDict, total=False):
        """Arguments for add_argument, with ``dest`` required."""

        action: _ActionStr | type[argparse.Action]
        nargs: int | _NArgsStr | _SUPPRESS_T | None
        const: Any
        default: Any
        type: Callable[[Any], Any] | argparse.FileType
        choices: Iterable[Any] | None
        required: bool
        help: str | None
        metavar: str | tuple[str, ...] | None
        dest: Required[str]
        version: str

    AnyArgumentArgs: TypeAlias = "AddArgumentArgs"


class PluginBase:
    """Optional base class for snakebids plugins providing utility methods."""

    PREFIX = ""

    @overload
    def pop(self, mapping: dict[str, Any], key: str, default: Any, /) -> Any: ...

    @overload
    def pop(self, mapping: dict[str, Any], key: str, /) -> Any: ...

    def pop(self, mapping: dict[str, Any], *args: Any):
        """Remove specified key from mapping, prepending the plugin prefix."""
        if not args:
            msg = "pop expected at least two arguments, got 1"
            raise TypeError(msg)
        if self.PREFIX:
            return mapping.pop(f"{self.PREFIX}.{args[0]}", *args[1:])
        return mapping.pop(*args)

    def get(self, mapping: dict[str, Any], key: str):
        """Retrieve specified key from mapping, prepending the plugin prefix."""
        if self.PREFIX:
            return mapping.get(f"{self.PREFIX}.{key}")
        return mapping.get(key)

    def assign(self, mapping: dict[str, Any], key: str, value: Any):
        """Assign specified key from mapping, prepending the plugin prefix."""
        if self.PREFIX:
            mapping[f"{self.PREFIX}.{key}"] = value
        else:
            mapping[key] = value

    def try_add_argument(
        self,
        parser: argparse.ArgumentParser | _ArgumentGroup,
        *name_or_flags: str,
        **kwargs: Unpack[AnyArgumentArgs],
    ) -> argparse.Action | None:
        """Add argument to parser if provided dest is not already present."""
        if self.PREFIX and not kwargs["dest"].startswith(self.PREFIX):
            kwargs["dest"] = f"{self.PREFIX}.{kwargs['dest']}"
        if kwargs["dest"] in {a.dest for a in parser._actions}:  # noqa: SLF001
            return None
        return parser.add_argument(*name_or_flags, **kwargs)

    def add_argument(
        self,
        parser: argparse.ArgumentParser | _ArgumentGroup,
        *name_or_flags: str,
        **kwargs: Unpack[AnyArgumentArgs],
    ) -> argparse.Action | None:
        """Add argument to parser, applying prefix to dest as necesary."""
        if self.PREFIX and not kwargs["dest"].startswith(self.PREFIX):
            kwargs["dest"] = f"{self.PREFIX}.{kwargs['dest']}"
        return parser.add_argument(*name_or_flags, **kwargs)
