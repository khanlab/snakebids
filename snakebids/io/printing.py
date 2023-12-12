from __future__ import annotations

import itertools as it
import json
import operator as op
import textwrap
from math import ceil, floor, inf
from typing import Sequence

import more_itertools as itx

from snakebids.types import ZipList


def quote_wrap(val: str) -> str:
    """Wrap string in quotes, with additional character escapes for printing."""
    return (
        json.dumps(val, ensure_ascii=False)
        .replace("\x85", "\\x85")
        .replace("\u2028", "\\u2028")
        .replace("\u2029", "\\u2029")
    )


def format_zip_lists(
    zip_list: ZipList, max_width: int | float | None = None, tabstop: int = 4
) -> str:
    """Pretty-format zip-lists in a tablar format.

    Parameters
    ----------
    zip_list
        zip_list to format
    max_width
        Maximum character width. If possible, zip_list values will be elided to fit
        within this width
    tabstop
        Number of spaces to include in each level of indentation
    """
    if not zip_list:
        return "{}"
    table = [_format_zip_row(key, row) for key, row in zip_list.items()]
    widths = [max(len(val) for val in col) for col in zip(*table)]
    aligned = _align_zip_table(table, widths)
    cols: list[Sequence[str]] = list(zip(*aligned))
    keys = cols[:2]
    vals = cols[2:]
    elided_cols = it.chain(
        keys,
        _elide_zip_table(
            vals,
            len(aligned),
            widths[2:],
            # max width reduces for the indent, minus 2 for the closing "]," minus the
            # width of the key and opening brace
            max_width=(max_width or inf) - tabstop - 2 - sum(widths[:2]),
        ),
        [["],\n"] * len(table)],
    )
    return "".join(
        [
            "{\n",
            textwrap.indent("".join(itx.flatten(zip(*elided_cols))), " " * tabstop),
            "}",
        ]
    )


def _format_zip_row(key: str, row: list[str]) -> list[str]:
    formatted_values = [
        quote_wrap(val) + sep
        for val, sep in it.zip_longest(row, it.repeat(", ", len(row) - 1), fillvalue="")
    ]
    row = [f"{quote_wrap(key)}: "]
    row.append("[")
    row.extend(formatted_values)
    return row


def _align_zip_table(table: list[list[str]], widths: list[int]) -> list[list[str]]:
    output: list[list[str]] = []
    for row in table:
        spaces = [" " * (width - len(val)) for val, width in zip(row, widths)]
        output.append(list(it.starmap(op.add, zip(row, spaces))))
    return output


def _elide_zip_table(
    cols: list[Sequence[str]], col_len: int, widths: list[int], max_width: int | float
) -> list[Sequence[str]]:
    def new_col(val: str):
        return [[val] * col_len]

    overflow = int(max(sum(widths) - max_width, 0))
    elision = _find_elision(list(widths), slice(0, 0), overflow)
    if elision != slice(0, 0):
        return list(
            it.chain(
                cols[: elision.start],
                new_col("..."),
                (new_col(" ") if elision.stop - elision.start < len(cols) else []),
                cols[elision.stop :],
            )
        )
    return cols


def _find_elision(widths: list[int], excluded: slice, overflow: int) -> slice:
    # Add 4 to overflow to account for elipses
    if max(sum(widths[excluded]) - 4, 0) >= overflow:
        return excluded
    span = excluded.stop - excluded.start
    # don't let span become longer than list of values
    if span >= len(widths):
        return excluded

    num_vals = len(widths)
    mid = floor(num_vals / 2)

    # need different rules for handling exclusions of even length depending on whether
    # theres an even or odd total number of values.
    left_bias = floor if num_vals % 2 else ceil
    right_bias = ceil if num_vals % 2 else floor

    return _find_elision(
        widths,
        slice(mid - left_bias(span / 2), mid + right_bias(span / 2) + 1, 1),
        overflow,
    )
