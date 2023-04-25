from __future__ import annotations

import pytest

from snakebids.core.filtering import filter_list
from snakebids.types import ZipLists


@pytest.mark.parametrize(
    ("zip_list", "filters", "output"),
    [
        (
            {
                "dir": ["AP", "PA", "AP", "PA", "AP", "PA", "AP", "PA"],
                "acq": ["98", "98", "98", "98", "99", "99", "99", "99"],
                "subject": ["01", "01", "02", "02", "01", "01", "02", "02"],
            },
            {"subject": "01"},
            {
                "dir": ["AP", "PA", "AP", "PA"],
                "acq": ["98", "98", "99", "99"],
                "subject": ["01", "01", "01", "01"],
            },
        ),
        (
            {
                "dir": ["AP", "PA", "AP", "PA", "AP", "PA", "AP", "PA"],
                "acq": ["98", "98", "98", "98", "99", "99", "99", "99"],
                "subject": ["01", "01", "02", "02", "01", "01", "02", "02"],
            },
            {"acq": "98"},
            {
                "dir": ["AP", "PA", "AP", "PA"],
                "acq": ["98", "98", "98", "98"],
                "subject": ["01", "01", "02", "02"],
            },
        ),
        (
            {
                "dir": ["AP", "PA", "AP", "PA", "AP", "PA", "AP", "PA"],
                "acq": ["98", "98", "98", "98", "99", "99", "99", "99"],
                "subject": ["01", "01", "02", "02", "01", "01", "02", "02"],
            },
            {"dir": "AP"},
            {
                "dir": ["AP", "AP", "AP", "AP"],
                "acq": ["98", "98", "99", "99"],
                "subject": ["01", "02", "01", "02"],
            },
        ),
        (
            {
                "dir": ["AP", "PA", "AP", "PA", "AP", "PA", "AP", "PA"],
                "acq": ["98", "98", "98", "98", "99", "99", "99", "99"],
                "subject": ["01", "01", "02", "02", "01", "01", "02", "02"],
            },
            {"subject": "03", "dir": "AP"},
            {"dir": [], "acq": [], "subject": []},
        ),
        (
            {
                "dir": ["AP", "PA", "AP", "PA", "AP", "PA", "AP", "PA"],
                "subject": ["01", "01", "02", "02", "03", "03", "04", "04"],
            },
            {"subject": ["01", "03"]},
            {
                "dir": ["AP", "PA", "AP", "PA"],
                "subject": ["01", "01", "03", "03"],
            },
        ),
        (
            {
                "dir": ["AP", "PA", "AP", "PA", "AP", "PA", "AP", "PA"],
                "subject": ["01", "01", "02", "02", "03", "03", "04", "04"],
            },
            {"subject": ["01", "02"], "dir": "AP"},
            {
                "dir": [
                    "AP",
                    "AP",
                ],
                "subject": ["01", "02"],
            },
        ),
    ],
)
def test_filter_list(
    zip_list: ZipLists,
    filters: dict[str, list[str] | str],
    output: dict[str, dict[str, list[str]]],
):
    assert filter_list(zip_list, filters) == output
