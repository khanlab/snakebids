from __future__ import absolute_import

from pathlib import Path

from snakebids.core.construct_bids import bids


def test_bids_subj():
    assert bids(root="bids", subject="001", suffix="T1w.nii.gz") == (
        "bids/sub-001/sub-001_T1w.nii.gz"
    )
    assert bids(root=Path("bids").resolve(), subject="001", suffix="T1w.nii.gz") == (
        str(Path.cwd() / "bids/sub-001/sub-001_T1w.nii.gz")
    )
