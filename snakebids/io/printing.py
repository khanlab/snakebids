from __future__ import annotations

import itertools as it
import json
import operator as op
import textwrap
from math import ceil, floor, inf

import more_itertools as itx

from snakebids.types import ZipList


def quote_wrap(val: str) -> str:
    return json.dumps(val, ensure_ascii=False)


def format_zip_lists(
    zip_list: ZipList, max_width: int | float | None = None, tabstop: int = 4
) -> str:
    table = [_format_zip_row(key, row) for key, row in zip_list.items()]
    widths = [max(len(val) for val in col) for col in zip(*table)]
    aligned = _align_zip_table(table, widths)
    elided_cols = it.chain(
        # max width reduces for the indent, plus 2 for the closing "]," on each line
        _elide_zip_table(aligned, widths, max_width=(max_width or inf) - tabstop - 2),
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
    return [f"{quote_wrap(key)}: "] + ["[" + formatted_values[0]] + formatted_values[1:]


def _align_zip_table(table: list[list[str]], widths: list[int]) -> list[list[str]]:
    output: list[list[str]] = []
    for row in table:
        spaces = [" " * (width - len(val)) for val, width in zip(row, widths)]
        output.append(list(it.starmap(op.add, zip(row, spaces))))
    return output


def _elide_zip_table(
    table: list[list[str]], widths: list[int], max_width: int | float
) -> list[list[str]] | list[tuple[str]] | list[list[str] | tuple[str]]:
    def new_col(val: str):
        return [[val] * len(table)]

    overflow = int(max(sum(widths) - (max_width or inf), 0))
    elision = _find_elision(list(widths), slice(0, 0), overflow)
    cols = list(zip(*table))
    if elision != slice(0, 0):
        return list(
            it.chain(
                cols[: elision.start]
                if elision.start > 1
                else [cols[0], *new_col("[")],
                new_col("..."),
                (new_col(" ") if elision.stop - elision.start < len(cols) - 1 else []),
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
    if span >= len(widths) - 1:
        return excluded

    # subtract 1 to account for the key, then add it back when getting the mid index
    num_vals = len(widths) - 1
    mid = floor(num_vals / 2) + 1

    # need different rules for handling exclusions of even length depending on whether
    # theres an even or odd total number of values.
    left_bias = floor if num_vals % 2 else ceil
    right_bias = ceil if num_vals % 2 else floor

    return _find_elision(
        widths,
        slice(mid - left_bias(span / 2), mid + right_bias(span / 2) + 1, 1),
        overflow,
    )
