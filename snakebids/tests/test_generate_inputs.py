import os

from bids import BIDSLayout
import pytest

from .. import generate_inputs, __get_lists_from_bids


def test_t1w():
    # create config
    real_bids_dir = "snakebids/tests/data/bids_t1w"
    derivatives = False
    pybids_inputs = {
        "t1": {
            "filters": {"suffix": "T1w"},
            "wildcards": ["acquisition", "subject", "session", "run"],
        }
    }

    # Can't define particpant_label and exclude_participant_label
    with pytest.raises(ValueError) as v_error:
        config = generate_inputs(
            pybids_inputs=pybids_inputs,
            bids_dir=real_bids_dir,
            derivatives=derivatives,
            participant_label="001",
            exclude_participant_label="002",
        )
        assert v_error.msg == (
            "Cannot define both participant_label and "
            "exclude_participant_label at the same time"
        )

    # Simplest case -- one input type, using pybids
    config = generate_inputs(
        pybids_inputs=pybids_inputs,
        bids_dir=real_bids_dir,
        derivatives=derivatives,
    )
    # Order of the subjects is not deterministic
    assert config["input_lists"] in [
        {"t1": {"acq": ["mprage"], "subject": ["001", "002"]}},
        {"t1": {"acq": ["mprage"], "subject": ["002", "001"]}},
    ]
    assert config["input_zip_lists"] == {
        "t1": {"acq": ["mprage", "mprage"], "subject": ["001", "002"]}
    }
    assert config["input_wildcards"] == {
        "t1": {"acq": "{acq}", "subject": "{subject}"}
    }
    assert config["subjects"] == ["001", "002"]
    assert config["sessions"] == []
    assert config["subj_wildcards"] == {"subject": "{subject}"}

    pybids_inputs_suffix = {
        "scan": {
            "filters": {},
            "wildcards": [
                "acquisition",
                "subject",
                "session",
                "run",
                "suffix",
            ],
        }
    }
    config = generate_inputs(
        pybids_inputs=pybids_inputs_suffix,
        bids_dir=real_bids_dir,
        derivatives=derivatives,
        participant_label="001",
    )
    assert config["input_lists"] == {
        "scan": {"acq": ["mprage"], "subject": ["001"], "suffix": ["T1w"]}
    }
    assert config["input_zip_lists"] == {
        "scan": {"acq": ["mprage"], "subject": ["001"], "suffix": ["T1w"]}
    }
    assert config["input_wildcards"] == {
        "scan": {"acq": "{acq}", "subject": "{subject}", "suffix": "{suffix}"}
    }
    assert config["subjects"] == ["001"]
    assert config["sessions"] == []
    assert config["subj_wildcards"] == {"subject": "{subject}"}

    # Two input types, specified by pybids or path override
    wildcard_path_t1 = os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        "data/bids_t1w",
        "sub-{subject}/anat/sub-{subject}_acq-{acq}_T1w.nii.gz",
    )
    wildcard_path_t2 = os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        "data/bids_t1w",
        "sub-{subject}/anat/sub-{subject}_T2w.nii.gz",
    )
    pybids_inputs = {
        "t1": {
            "filters": {"suffix": "T1w"},
            "wildcards": ["acquisition", "subject", "session", "run"],
        },
        "t2": {
            "filters": {"suffix": "T2w"},
            "wildcards": ["acquisition", "subject", "session", "run"],
        },
    }
    bids_dir = real_bids_dir
    for idx in range(2):
        if idx == 1:
            pybids_inputs["t1"]["custom_path"] = wildcard_path_t1
        if idx == 2:
            pybids_inputs["t2"]["custom_path"] = wildcard_path_t2
            bids_dir = "-"
        config = generate_inputs(
            pybids_inputs=pybids_inputs,
            bids_dir=bids_dir,
            derivatives=derivatives,
        )
        # Order of the subjects is not deterministic
        assert config["input_lists"]["t1"] in [
            {"acq": ["mprage"], "subject": ["001", "002"]},
            {"acq": ["mprage"], "subject": ["002", "001"]},
        ]
        assert config["input_lists"]["t2"] == {"subject": ["002"]}
        assert config["input_zip_lists"]["t1"] in [
            {"acq": ["mprage", "mprage"], "subject": ["001", "002"]},
            {"acq": ["mprage", "mprage"], "subject": ["002", "001"]},
        ]
        assert config["input_zip_lists"]["t2"] == {"subject": ["002"]}
        assert config["input_wildcards"] == {
            "t1": {"acq": "{acq}", "subject": "{subject}"},
            "t2": {"subject": "{subject}"},
        }
        assert config["subjects"] == ["001", "002"]
        assert config["sessions"] == []
        assert config["subj_wildcards"] == {"subject": "{subject}"}


def test_get_lists_from_bids():
    bids_dir = "snakebids/tests/data/bids_t1w"
    wildcard_path_t1 = os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        "data/bids_t1w",
        "sub-{subject}/anat/sub-{subject}_acq-{acq}_T1w.nii.gz",
    )
    wildcard_path_t2 = os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        "data/bids_t1w",
        "sub-{subject}/anat/sub-{subject}_T2w.nii.gz",
    )
    print(wildcard_path_t1)
    layout = BIDSLayout(bids_dir, validate=False)
    pybids_inputs = {
        "t1": {
            "filters": {"suffix": "T1w"},
            "wildcards": ["acquisition", "subject", "session", "run"],
        },
        "t2": {
            "filters": {"suffix": "T2w"},
            "wildcards": ["acquisition", "subject", "session", "run"],
        },
    }

    # Want to test both inputs from layout, both inputs from custom path, and
    # one of each. This setup should produce the same results every time.
    for idx in range(2):
        if idx == 1:
            pybids_inputs["t1"]["custom_path"] = wildcard_path_t1
        elif idx == 2:
            pybids_inputs["t2"]["custom_path"] = wildcard_path_t2

        config = __get_lists_from_bids(layout, pybids_inputs)
        assert config["input_path"] == {
            "t1": wildcard_path_t1,
            "t2": wildcard_path_t2,
        }
        assert config["input_zip_lists"]["t1"] in [
            {"acq": ["mprage", "mprage"], "subject": ["001", "002"]},
            {"acq": ["mprage", "mprage"], "subject": ["002", "001"]},
        ]
        assert config["input_zip_lists"]["t2"] == {"subject": ["002"]}
        # The order of multiple wildcard values is not deterministic
        assert config["input_lists"] in [
            {
                "t1": {"acq": ["mprage"], "subject": ["001", "002"]},
                "t2": {"subject": ["002"]},
            },
            {
                "t1": {"acq": ["mprage"], "subject": ["002", "001"]},
                "t2": {"subject": ["002"]},
            },
        ]
        assert config["input_wildcards"] == {
            "t1": {"acq": "{acq}", "subject": "{subject}"},
            "t2": {"subject": "{subject}"},
        }
