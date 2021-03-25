"""Tests for snakemake_io"""

import collections

from .. import snakemake_io


def test_glob_wildcards():
    """Test glob_wildcards() with various patterns"""
    file_path = (
        "snakebids/tests/data/bids_t1w/sub-001/anat/"
        "sub-001_acq-mprage_T1w.nii.gz"
    )
    empty_wildcards = collections.namedtuple("Wildcards", [])()
    assert snakemake_io.glob_wildcards(file_path) == empty_wildcards

    acq_wildcard_path = (
        "snakebids/tests/data/bids_t1w/sub-001/anat/"
        "sub-001_acq-{acq}_T1w.nii.gz"
    )
    acq_wildcards = collections.namedtuple("Wildcards", ["acq"])(["mprage"])
    assert snakemake_io.glob_wildcards(acq_wildcard_path) == acq_wildcards

    # Order of wildcards in the namedtuple is not deterministic
    both_wildcard_path = (
        "snakebids/tests/data/bids_t1w/sub-{subject}/anat/"
        "sub-{subject}_acq-{acq}_T1w.nii.gz"
    )
    both_wildcards = [
        collections.namedtuple("Wildcards", ["subject", "acq"])(
            ["001", "002"], ["mprage", "mprage"]
        ),
        collections.namedtuple("Wildcards", ["subject", "acq"])(
            ["002", "001"], ["mprage", "mprage"]
        ),
        collections.namedtuple("Wildcards", ["acq", "subject"])(
            ["mprage", "mprage"], ["001", "002"]
        ),
        collections.namedtuple("Wildcards", ["acq", "subject"])(
            ["mprage", "mprage"], ["002", "001"]
        ),
    ]
    both_wildcards_one_file = [
        collections.namedtuple("Wildcards", ["subject", "acq"])(
            ["001"], ["mprage"]
        ),
        collections.namedtuple("Wildcards", ["acq", "subject"])(
            ["mprage"], ["001"]
        ),
    ]

    assert snakemake_io.glob_wildcards(both_wildcard_path) in both_wildcards
    assert (
        snakemake_io.glob_wildcards(both_wildcard_path, files=[file_path])
        in both_wildcards_one_file
    )
