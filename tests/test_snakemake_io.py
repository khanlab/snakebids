"""Tests for snakemake_io"""

from snakebids.utils import snakemake_io


def test_glob_wildcards():
    """Test glob_wildcards() with various patterns"""
    file_path = "tests/data/bids_t1w/sub-001/anat/sub-001_acq-mprage_T1w.nii.gz"
    empty_wildcards = {}
    assert snakemake_io.glob_wildcards(file_path) == empty_wildcards

    acq_wildcard_path = "tests/data/bids_t1w/sub-001/anat/sub-001_acq-{acq}_T1w.nii.gz"
    acq_wildcards = {"acq": ["mprage"]}
    assert snakemake_io.glob_wildcards(acq_wildcard_path) == acq_wildcards

    # Order of wildcards in the namedtuple is not deterministic
    both_wildcard_path = (
        "tests/data/bids_t1w/sub-{subject}/anat/sub-{subject}_acq-{acq}_T1w.nii.gz"
    )
    both_wildcards = [
        {"subject": ["001", "002"], "acq": ["mprage", "mprage"]},
        {"subject": ["002", "001"], "acq": ["mprage", "mprage"]},
        {"acq": ["mprage", "mprage"], "subject": ["001", "002"]},
        {"acq": ["mprage", "mprage"], "subject": ["002", "001"]},
    ]
    both_wildcards_one_file = [
        {"subject": ["001"], "acq": ["mprage"]},
        {"acq": ["mprage"], "subject": ["001"]},
    ]

    assert snakemake_io.glob_wildcards(both_wildcard_path) in both_wildcards
    assert (
        snakemake_io.glob_wildcards(both_wildcard_path, files=[file_path])
        in both_wildcards_one_file
    )


def test_glob_wildcards_includes_unmatched_wildcards():
    assert snakemake_io.glob_wildcards("tpl-{tpl}/{suffix}.nii", files=[]) == {
        "tpl": [],
        "suffix": [],
    }


def test_glob_wildcards_normalizes_files():
    assert snakemake_io.glob_wildcards(
        "https://example.com/tpl-{tpl}/{suffix}.nii",
        files=["https://example.com/tpl-A/T1w.nii", "https://example.com/tpl-A/b/c"],
    ) == {"tpl": ["A"], "suffix": ["T1w"]}
    assert snakemake_io.glob_wildcards(
        "tpl-{tpl}/{suffix}.nii", files=["./tpl-A/T1w.nii"]
    ) == {"tpl": ["A"], "suffix": ["T1w"]}
