from __future__ import annotations

import itertools as it
from collections.abc import Iterable
from typing import Any

import more_itertools as itx

from snakebids.types import ZipListLike
from snakebids.utils.snakemake_templates import MissingEntityError, SnakemakeFormatter

try:
    from snakemake.io import AnnotatedString as _AnnotatedString  # type: ignore

except ImportError:

    def _get_flags(path: Any) -> dict[str, Any] | None:
        """Return None when snakemake is not installed (no AnnotatedString support)."""
        return None

    def _make_annotated(s: str, flags: dict[str, Any]) -> str:
        """Return the string unchanged when snakemake is not installed."""
        return s
else:

    def _get_flags(path: Any) -> dict[str, Any] | None:
        """Extract flags from an AnnotatedString path."""
        if isinstance(path, _AnnotatedString) and path.flags:  # type: ignore
            return dict(path.flags)  # type: ignore
        return None

    def _make_annotated(s: str, flags: dict[str, Any]) -> str:
        """Create an AnnotatedString with the given flags copied from the template."""
        result = _AnnotatedString(s)
        result.flags.update(flags)  # type: ignore
        return result


def expand(
    paths: Iterable[Any],
    zip_lists: ZipListLike,
    allow_missing: bool,
    **wildcards: Any,
) -> list[str]:
    """Expand template paths using SnakemakeFormatter.

    For each template path, iterates over rows of the zip-list and optional
    extra wildcard combinations (product), formatting each path per row.
    None values in wildcards are converted to "".
    Output is order-preserving deduplicated.
    """
    formatter = SnakemakeFormatter(allow_missing=allow_missing)

    # Normalize extra wildcards: convert None (scalar or in list) → "", ensure lists
    extra: dict[str, list[str]] = {}
    for k, vals in wildcards.items():
        if vals is None:
            extra[k] = [""]
        else:
            extra[k] = ["" if v is None else v for v in itx.always_iterable(vals)]

    rows = list(zip(*zip_lists.values(), strict=True))

    results: list[str] = []
    for path in paths:
        path_str = str(path)
        flags = _get_flags(path)
        for row, *combo in it.product(rows, *extra.values()):
            kwargs = {
                **dict(zip(zip_lists, row, strict=True)),
                **dict(zip(extra, combo, strict=True)),
            }
            try:
                result = formatter.vformat(path_str, (), kwargs)
            except MissingEntityError as err:
                msg = f"no values given for wildcard {err.entity!r}."
                raise KeyError(msg) from err
            except KeyError as err:
                msg = f"no values given for wildcard {err.args[0]!r}."
                raise KeyError(msg) from err

            if flags is not None:
                result = _make_annotated(result, flags)
            results.append(result)

    return list(dict.fromkeys(results))
