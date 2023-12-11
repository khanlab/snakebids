from __future__ import annotations

import itertools as it
import os
import sys
from pathlib import Path
from typing import Protocol

import more_itertools as itx

from snakebids.paths._utils import BidsPathSpec, find_entity


class BidsFunction(Protocol):
    """Signature for functions returned by ``bids_factory``.

    See :func:`bids` for more details
    """

    def __call__(
        self,
        root: str | Path | None = None,
        *,
        datatype: str | None = None,
        prefix: str | None = None,
        suffix: str | None = None,
        extension: str | None = None,
        **entities: str | bool,
    ) -> str:
        ...


def _get_entity_parser(aliases: dict[str, str]):
    def parse_entities(entities: dict[str, str | bool]) -> dict[str, str]:
        result: dict[str, str] = {}
        for entity, val in entities.items():
            # strip underscores from keys (needed so that users can use reserved
            # keywords by appending a _)
            stripped = entity.rstrip("_")
            unaliased = aliases.get(stripped, stripped)
            if unaliased in result:
                aliased = itx.nth(aliases, list(aliases.values()).index(unaliased))
                err = (
                    "Long and short names of an entity cannot be used in the same "
                    f"call to bids(): got '{aliased}' and '{unaliased}'"
                )
                raise ValueError(err)
            result[unaliased] = str(val)
        return result

    return parse_entities


def bids_factory(spec: BidsPathSpec, *, _v0_0_0: bool = False) -> BidsFunction:
    """Generate bids functions according to the supplied spec.

    Parameters
    ----------
        spec
            Valid Bids Spec object
        _v0_0_0
            Provides backward compatibility for the bids_v0_0_0 signature. Should not
            otherwise be used
    """
    order: list[str] = []
    dirs: set[str] = set()
    aliases: dict[str, str] = {}

    if _v0_0_0:
        subject_dir_default = find_entity(spec, "subject").get("dir", False)
        session_dir_default = find_entity(spec, "session").get("dir", False)
    else:
        subject_dir_default = True
        session_dir_default = True
    for entry in spec:
        tag = entry.get("tag", entry["entity"])
        order.append(tag)
        aliases[entry["entity"]] = tag
        if entry.get("dir"):
            dirs.add(tag)

    parse_entities = _get_entity_parser(aliases)

    def bids(
        root: str | Path | None = None,
        *,
        datatype: str | None = None,
        prefix: str | None = None,
        suffix: str | None = None,
        extension: str | None = None,
        **entities: str | bool,
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
        if _v0_0_0:
            from snakebids.paths.specs import v0_0_0

            include_subject_dir = bool(entities.pop("include_subject_dir", True))
            include_session_dir = bool(entities.pop("include_session_dir", True))
            if (
                include_session_dir ^ session_dir_default
                or include_subject_dir ^ subject_dir_default
            ):
                return bids_factory(v0_0_0(include_subject_dir, include_session_dir))(
                    root,
                    datatype=datatype,
                    prefix=prefix,
                    suffix=suffix,
                    extension=extension,
                    **entities,
                )

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

        parsed = parse_entities(entities)

        spec_parts: list[str] = []
        custom_parts: list[str] = []
        split: int = sys.maxsize + 1
        path_parts: list[str] = []

        if root:
            path_parts.append(str(root))
        if prefix:
            spec_parts.append(prefix)
        for entity in order:
            # Check for `*` first so that if user specifies an entity called `*` we
            # don't skip setting the split
            if entity == "*":
                split = len(spec_parts)
            elif value := parsed.pop(entity, None):
                spec_parts.append(f"{entity}-{value}")
                if entity in dirs:
                    path_parts.append(f"{entity}-{value}")
        for key, value in parsed.items():
            custom_parts.append(f"{key}-{value}")

        if datatype:
            path_parts.append(datatype)
        path_parts.append(
            "_".join(it.chain(spec_parts[:split], custom_parts, spec_parts[split:]))
        )

        tail = f"_{suffix}{extension or ''}" if suffix else extension or ""

        return os.path.join(*path_parts) + tail

    return bids
