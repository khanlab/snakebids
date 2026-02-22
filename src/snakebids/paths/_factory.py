from __future__ import annotations

import itertools as it
import os
import sys
import textwrap
import warnings
from collections.abc import Iterator
from pathlib import Path
from typing import Literal, Protocol, TypeAlias

import more_itertools as itx

from snakebids.io.console import in_interactive_session
from snakebids.paths import specs
from snakebids.paths._utils import (
    OPTIONAL_WILDCARD,
    BidsFlags,
    BidsPathSpec,
    find_entity,
)
from snakebids.utils.snakemake_templates import SnakemakeWildcards


class BidsFunction(Protocol):
    """Signature for functions returned by ``bids_factory``.

    See :func:`bids` for more details
    """

    def __call__(
        self,
        root: str | Path | None = None,
        *,
        datatype: str | BidsFlags | None = None,
        prefix: str | BidsFlags | None = None,
        suffix: str | BidsFlags | None = None,
        extension: str | BidsFlags | None = None,
        **entities: str | BidsFlags,
    ) -> str: ...


_WildcardKinds: TypeAlias = Literal["fixed", "optional"]


class _LeadingUnderscore:
    def __init__(self, kind: _WildcardKinds, emitter: _UnderscoreEmitter):
        self.emitter = emitter
        self.kind: _WildcardKinds = kind

    def __str__(self) -> str:
        return self.emitter.next(self.kind)


class _UnderscoreEmitter:
    def __init__(self) -> None:
        self.emitter = None
        self.emitter_kinds: dict[_WildcardKinds, Iterator[str]] = {
            "fixed": it.repeat("_"),
            "optional": it.chain([str(SnakemakeWildcards.underscore)], it.repeat("_")),
        }

    def new(self, kind: _WildcardKinds):
        return _LeadingUnderscore(kind, self)

    def next(self, kind: _WildcardKinds):
        if self.emitter is None:
            self.emitter = self.emitter_kinds[kind]
            return ""
        if kind == "optional":
            return ""
        return next(self.emitter)


class _NameBuilder(list[str | _LeadingUnderscore]):
    def __init__(self, emitter: _UnderscoreEmitter):
        self.emitter = emitter

    def append_entity(self, key: str, value: str | BidsFlags) -> None:
        if value is OPTIONAL_WILDCARD:
            wcard_kind = "optional"
            wc = SnakemakeWildcards(key)
            value = f"{wc.dummy}{wc.variable}"
        else:
            wcard_kind = "fixed"
            value = f"{key}-{value}"

        self.append(self.emitter.new(wcard_kind))
        self.append(value)

    def append_literal(self, literal: str) -> None:
        self.append(self.emitter.new("fixed"))
        self.append(literal)

    def append_suffix(self, suffix: str | BidsFlags) -> None:
        self.append_literal(
            str(SnakemakeWildcards.suffix) if suffix is OPTIONAL_WILDCARD else suffix
        )

    def append_extension(self, extension: str | BidsFlags) -> None:
        self.append(
            str(SnakemakeWildcards.extension)
            if extension is OPTIONAL_WILDCARD
            else extension
        )

    def sanitized(self):
        return list(filter(lambda s: isinstance(s, str), self))


class _PathBuilder(list[str]):
    def append_entity(self, key: str, value: str | BidsFlags):
        self.append(
            str(SnakemakeWildcards(key).directory)
            if value is OPTIONAL_WILDCARD
            else f"{key}-{value}{os.sep}"
        )

    def append_datatype(self, datatype: str | BidsFlags):
        self.append(
            f"{SnakemakeWildcards.datatype}{SnakemakeWildcards.slash}"
            if datatype is OPTIONAL_WILDCARD
            else datatype + os.sep
        )


def bids_factory(spec: BidsPathSpec, *, _implicit: bool = False) -> BidsFunction:
    """Generate bids functions according to the supplied spec.

    Parameters
    ----------
        spec
            Valid Bids Spec object
        _implicit
            Flag used internally to mark the default generated bids function. The
            resulting builder will warn when custom entities are used
    """
    return _Bids(spec, _implicit=_implicit)


class _Bids:
    def __init__(self, spec: BidsPathSpec, *, _implicit: bool = False):
        self.spec = spec
        self._implicit = _implicit

        self.order: list[str] = []
        self.dirs: set[str] = set()
        self.aliases: dict[str, str] = {}

        self.subject_dir_default = find_entity(spec, "subject").get("dir", False)
        self.session_dir_default = find_entity(spec, "session").get("dir", False)

        for entry in spec:
            tag = entry.get("tag", entry["entity"])
            self.order.append(tag)
            self.aliases[entry["entity"]] = tag
            if entry.get("dir"):
                self.dirs.add(tag)

    def parse_entities(self, entities: dict[str, str | BidsFlags]):
        result: dict[str, str | BidsFlags] = {}
        for entity, val in entities.items():
            # strip underscores from keys (needed so that users can use reserved
            # keywords by appending a _)
            stripped = entity.rstrip("_")
            unaliased = self.aliases.get(stripped, stripped)
            if unaliased in result:
                aliased = itx.nth(
                    self.aliases, list(self.aliases.values()).index(unaliased)
                )
                err = (
                    "Long and short names of an entity cannot be used in the same "
                    f"call to bids(): got '{aliased}' and '{unaliased}'"
                )
                raise ValueError(err)
            # Check if the value is the OPTIONAL_WILDCARD enum
            result[unaliased] = val
        return result

    def handle_subses_dir(
        self,
        root: str | Path | None,
        /,
        *,
        sub_dir: str | bool | BidsFlags | None,
        ses_dir: str | bool | BidsFlags | None,
        datatype: str | BidsFlags | None = None,
        prefix: str | BidsFlags | None = None,
        suffix: str | BidsFlags | None = None,
        extension: str | BidsFlags | None = None,
        **entities: str | BidsFlags,
    ) -> str | None:
        newspec = None
        warn_subses_dir = False

        if isinstance(ses_dir, bool):
            warn_subses_dir = True
            if ses_dir ^ self.session_dir_default:
                if newspec is None:
                    newspec = self.spec.copy()
                find_entity(newspec, "session")["dir"] = ses_dir

        if isinstance(sub_dir, bool):
            warn_subses_dir = True
            if sub_dir ^ self.subject_dir_default:
                if newspec is None:
                    newspec = self.spec.copy()
                find_entity(newspec, "subject")["dir"] = sub_dir

        if warn_subses_dir:
            wrn_msg = (
                "include_session_dir and include_subject_dir are deprecated and "
                "will be removed in a future release. Builder functions without "
                "directories can be created using the bids_factory and spec "
                "functions:\n"
                "\tfrom snakebids.paths import bids_factory, specs\n"
                "\tbids_ses = bids_factory(specs.v0_0_0(subject_dir=False))\n"
                "\tbids_ses(...)\n"
            ).expandtabs(4)
            warnings.warn(wrn_msg, stacklevel=4)
        if newspec:
            return _Bids(newspec)(
                root,
                datatype=datatype,
                prefix=prefix,
                suffix=suffix,
                extension=extension,
                **entities,
            )
        return None

    def check_custom_parts_without_explicit_spec(
        self, custom_parts: _NameBuilder, path: str
    ):
        if custom_parts and self._implicit and not in_interactive_session():
            wrn_msg = (
                textwrap.fill(
                    "Path generated with unrecognized entities, and a snakebids spec "
                    "has not been explicitly declared. This could break in future "
                    "snakebids versions, as the default spec can be changed without "
                    "warning.",
                ),
                f"\tpath = {path!r}",
                f"\tentities = {custom_parts.sanitized()!r}",
                "",
                "Please declare a spec using:",
                "\tfrom snakebids import set_bids_spec",
                f'\tset_bids_spec("{specs.LATEST}")',
            )
            warnings.warn("\n".join(wrn_msg), stacklevel=3)

    def __call__(  # noqa: PLR0912
        self,
        root: str | Path | None = None,
        *,
        datatype: str | BidsFlags | None = None,
        prefix: str | BidsFlags | None = None,
        suffix: str | BidsFlags | None = None,
        extension: str | BidsFlags | None = None,
        **entities: str | BidsFlags,
    ) -> str:
        """Generate bids or bids-like paths.

        Path will have the general form::

            [root]/[sub-{{subject}}]/[ses-{{session}}]/
            [prefix]_[sub-{{subject}}]_[ses-{{session}}]_[{{key}}-{{val}}_...]_[suffix]

        If no arguments are specified, an empty string will be returned.

        Parameters
        ----------
        root
            Root folder to include in the path (e.g. ``results``)
        datatype
            Folder to include after sub-/ses- (e.g. ``anat``, ``dwi`` )
        prefix
            String to prepend to the file name. Useful for injecting custom entities at
            the front of the filename, e.g. ``tpl-{{tpl}}``
        suffix
            Suffix plus, optionally, the extension (e.g. ``T1w.nii.gz``)
        extension
            bids extension, beginning with ``.`` (e.g. ``.nii.gz``).  Typically
            shouldn't be specified manually: extensions should be listed along with the
            suffix.
        entities
            bids entities as keyword arguments paired with values (e.g. ``space="T1w"``
            for ``space-T1w``)
        """
        if (
            result := self.handle_subses_dir(
                root,
                sub_dir=entities.pop("include_subject_dir", None),
                ses_dir=entities.pop("include_session_dir", None),
                datatype=datatype,
                prefix=prefix,
                suffix=suffix,
                extension=extension,
                **entities,
            )
        ) is not None:
            return result

        if prefix is OPTIONAL_WILDCARD:
            msg = "prefix may not be specified as optional"
            raise ValueError(msg)

        if not any([entities, suffix, extension]) and any([datatype, prefix]):
            raise ValueError(
                "At least one of suffix, extension, or an entity must be "
                "supplied.\n\tGot only: "
                + " and ".join(
                    filter(
                        None,
                        (
                            f"datatype='{datatype}'" if datatype else None,
                            f"prefix='{prefix}'" if prefix else None,
                        ),
                    )
                )
            )

        parsed = self.parse_entities(entities)

        emitter = _UnderscoreEmitter()
        spec_parts = _NameBuilder(emitter)
        custom_parts = _NameBuilder(emitter)
        split: int = sys.maxsize
        path_parts = _PathBuilder()

        if prefix is not None:
            spec_parts.append_literal(prefix)

        for entity in self.order:
            # Check for `*` first so that if user specifies an entity called `*` we
            # don't skip setting the split
            if entity == "*":
                split = len(spec_parts)
                continue

            if value := parsed.pop(entity, None):
                spec_parts.append_entity(entity, value)
                if entity in self.dirs:
                    path_parts.append_entity(entity, value)

        split = min(split, len(spec_parts))

        for key, value in parsed.items():
            custom_parts.append_entity(key, value)

        if datatype is not None:
            path_parts.append_datatype(datatype)

        if suffix is not None:
            spec_parts.append_suffix(suffix)

        if extension is not None:
            spec_parts.append_extension(extension)

        all_parts = it.chain(
            path_parts,
            spec_parts[:split],
            custom_parts,
            spec_parts[split:],
        )

        result = "".join(map(str, all_parts))

        if root is not None:
            result = os.path.join(root, result)

        self.check_custom_parts_without_explicit_spec(custom_parts, result)

        return result
