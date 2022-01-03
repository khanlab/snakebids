from __future__ import absolute_import

from pathlib import Path

from .. import bids


def test_bids_subj():
    assert bids(root="bids", subject="001", suffix="T1w.nii.gz") == Path(
        "bids/sub-001/sub-001_T1w.nii.gz"
    )
    assert bids(root=Path("bids").resolve(), subject="001", suffix="T1w.nii.gz") == (
        Path.cwd() / "bids/sub-001/sub-001_T1w.nii.gz"
    )
