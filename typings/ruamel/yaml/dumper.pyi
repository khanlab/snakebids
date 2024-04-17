"""
This type stub file was generated by pyright.
"""

from typing import Any, Optional

from ruamel.yaml.emitter import Emitter
from ruamel.yaml.representer import (
    BaseRepresenter,
    Representer,
    RoundTripRepresenter,
    SafeRepresenter,
)
from ruamel.yaml.resolver import BaseResolver, Resolver, VersionedResolver
from ruamel.yaml.serializer import Serializer

if False: ...
__all__ = ["BaseDumper", "SafeDumper", "Dumper", "RoundTripDumper"]

class BaseDumper(Emitter, Serializer, BaseRepresenter, BaseResolver):
    def __init__(
        self: Any,
        stream: StreamType,
        default_style: Any = ...,
        default_flow_style: Any = ...,
        canonical: Optional[bool] = ...,
        indent: Optional[int] = ...,
        width: Optional[int] = ...,
        allow_unicode: Optional[bool] = ...,
        line_break: Any = ...,
        encoding: Any = ...,
        explicit_start: Optional[bool] = ...,
        explicit_end: Optional[bool] = ...,
        version: Any = ...,
        tags: Any = ...,
        block_seq_indent: Any = ...,
        top_level_colon_align: Any = ...,
        prefix_colon: Any = ...,
    ) -> None: ...

class SafeDumper(Emitter, Serializer, SafeRepresenter, Resolver):
    def __init__(
        self,
        stream: StreamType,
        default_style: Any = ...,
        default_flow_style: Any = ...,
        canonical: Optional[bool] = ...,
        indent: Optional[int] = ...,
        width: Optional[int] = ...,
        allow_unicode: Optional[bool] = ...,
        line_break: Any = ...,
        encoding: Any = ...,
        explicit_start: Optional[bool] = ...,
        explicit_end: Optional[bool] = ...,
        version: Any = ...,
        tags: Any = ...,
        block_seq_indent: Any = ...,
        top_level_colon_align: Any = ...,
        prefix_colon: Any = ...,
    ) -> None: ...

class Dumper(Emitter, Serializer, Representer, Resolver):
    def __init__(
        self,
        stream: StreamType,
        default_style: Any = ...,
        default_flow_style: Any = ...,
        canonical: Optional[bool] = ...,
        indent: Optional[int] = ...,
        width: Optional[int] = ...,
        allow_unicode: Optional[bool] = ...,
        line_break: Any = ...,
        encoding: Any = ...,
        explicit_start: Optional[bool] = ...,
        explicit_end: Optional[bool] = ...,
        version: Any = ...,
        tags: Any = ...,
        block_seq_indent: Any = ...,
        top_level_colon_align: Any = ...,
        prefix_colon: Any = ...,
    ) -> None: ...

class RoundTripDumper(Emitter, Serializer, RoundTripRepresenter, VersionedResolver):
    def __init__(
        self,
        stream: StreamType,
        default_style: Any = ...,
        default_flow_style: Optional[bool] = ...,
        canonical: Optional[int] = ...,
        indent: Optional[int] = ...,
        width: Optional[int] = ...,
        allow_unicode: Optional[bool] = ...,
        line_break: Any = ...,
        encoding: Any = ...,
        explicit_start: Optional[bool] = ...,
        explicit_end: Optional[bool] = ...,
        version: Any = ...,
        tags: Any = ...,
        block_seq_indent: Any = ...,
        top_level_colon_align: Any = ...,
        prefix_colon: Any = ...,
    ) -> None: ...
