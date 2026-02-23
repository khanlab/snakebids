from __future__ import annotations

import itertools as it
from collections.abc import Iterable
from typing import Any

import more_itertools as itx

from snakebids.snakemake_compat import WildcardError
from snakebids.types import ZipList
from snakebids.utils.snakemake_templates import SnakemakeFormatter

try:
    from snakemake.io import AnnotatedString as _AnnotatedString

    def _get_flags(path: Any) -> dict[str, Any] | None:
        """Extract flags from an AnnotatedString path for preservation during expansion."""
        if isinstance(path, _AnnotatedString) and path.flags:
            return dict(path.flags)
        return None

    def _make_annotated(s: str, flags: dict[str, Any]) -> str:
        """Create an AnnotatedString with the given flags copied from the template."""
        result = _AnnotatedString(s)
        result.flags.update(flags)
        return result

except ImportError:

    def _get_flags(path: Any) -> dict[str, Any] | None:  # type: ignore[misc]
        """Return None when snakemake is not installed (no AnnotatedString support)."""
        return None

    def _make_annotated(s: str, flags: dict[str, Any]) -> str:  # type: ignore[misc]
        """Return the string unchanged when snakemake is not installed."""
        return s


def expand(
    paths: Iterable[Any],
    zip_lists: ZipList,
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

    # Build zip-list rows
    if zip_lists:
        rows: list[dict[str, str]] = [
            dict(zip(zip_lists.keys(), vals, strict=True))
            for vals in zip(*zip_lists.values(), strict=True)
        ]
    else:
        rows = [{}]

    # Compute extra wildcard combinations (product)
    if extra:
        extra_keys = list(extra.keys())
        extra_combos: list[tuple[str, ...]] = list(it.product(*extra.values()))
    else:
        extra_keys = []
        extra_combos = [()]

    try:
        results: list[str] = []
        for path in paths:
            path_str = str(path)
            flags = _get_flags(path)
            for row in rows:
                for combo in extra_combos:
                    kwargs = {**row, **dict(zip(extra_keys, combo, strict=True))}
                    result: str = formatter.format(path_str, **kwargs)
                    if flags is not None:
                        result = _make_annotated(result, flags)
                    results.append(result)
    except KeyError as e:
        raise WildcardError(f"No values given for wildcard {e}.") from e

    return list(dict.fromkeys(results))
