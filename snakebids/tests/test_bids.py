from pathlib import Path
from .. import bids


def test_bids_subj():
    assert (
        bids(root="bids", subject="001", suffix="T1w.nii.gz")
        == "bids/sub-001/sub-001_T1w.nii.gz"
    )
    assert (
        bids(root=Path("bids"), subject="001", suffix="T1w.nii.gz")
        == str(Path.cwd() / "bids/sub-001/sub-001_T1w.nii.gz")
    )
