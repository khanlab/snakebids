"""Hypothesis strategies for snakemake template tests."""

from __future__ import annotations

from collections.abc import Collection

import hypothesis.strategies as st


def _format_strings(
    *,
    base_exclude: set[str],
    extra_exclude: Collection[str] | None = None,
    min_size: int = 0,
):
    if extra_exclude is not None:
        exclude = base_exclude | set(extra_exclude)
    else:
        exclude = base_exclude

    return st.text(st.characters(exclude_characters=exclude), min_size=min_size)


def field_names(
    *, exclude_characters: Collection[str] | None = None, min_size: int = 0
):
    return _format_strings(
        base_exclude=set("{},"), extra_exclude=exclude_characters, min_size=min_size
    )


def safe_field_names(
    *, exclude_characters: Collection[str] | None = None, min_size: int = 0
):
    return _format_strings(
        base_exclude=set("{},_[."), extra_exclude=exclude_characters, min_size=min_size
    )


def constraints(*, exclude_characters: Collection[str] | None = None):
    return _format_strings(base_exclude=set("{}"), extra_exclude=exclude_characters)


def literals(*, exclude_characters: Collection[str] | None = None):
    return st.text(st.characters(exclude_characters=exclude_characters)).map(
        lambda s: s.replace("{", "{{").replace("}", "}}")
    )
