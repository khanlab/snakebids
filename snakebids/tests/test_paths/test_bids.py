from __future__ import annotations

from pathlib import Path

import pytest

from snakebids.paths.presets import bids


def test_bids_subj():
    assert bids(root="bids", subject="001", suffix="T1w.nii.gz") == (
        "bids/sub-001/sub-001_T1w.nii.gz"
    )
    assert bids(root=Path("bids").resolve(), subject="001", suffix="T1w.nii.gz") == (
        str(Path.cwd() / "bids/sub-001/sub-001_T1w.nii.gz")
    )


def test_bids_with_no_args_gives_empty_path():
    assert not bids()


@pytest.mark.parametrize(
    ("args"),
    [
        {"datatype": "foo"},
        {"prefix": "foo"},
        {"datatype": "foo", "prefix": "foo"},
    ],
)
def test_missing_essential_entities_gives_error(args: dict[str, str]):
    with pytest.raises(ValueError):
        bids(**args)
