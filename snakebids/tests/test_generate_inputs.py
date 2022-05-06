from __future__ import absolute_import

import filecmp
import itertools as it
import os
import shutil
from collections import defaultdict
from pathlib import Path

import pytest
from bids import BIDSLayout

from snakebids.core.construct_bids import bids
from snakebids.core.input_generation import (
    _gen_bids_layout,
    _get_lists_from_bids,
    generate_inputs,
)
from snakebids.tests.helpers import BidsListCompare


class TestAbsentConfigEntries:
    def get_entities(self, root):
        # Generate directory
        entities = {"subject": ["001", "002"], "acq": sorted(["foo", "bar"])}
        zip_list = defaultdict(list)
        for e in it.product(*entities.values()):
            d = dict(zip(entities.keys(), e))
            for key, val in d.items():
                zip_list[key].append(val)
            path = Path(bids(root, datatype="anat", suffix="T1w.nii.gz", **d))
            path.parent.mkdir(parents=True, exist_ok=True)
            path.touch()
        return entities, zip_list

    def test_missing_filters(self, tmpdir: Path):
        entities, zip_list = self.get_entities(tmpdir)

        # create config
        derivatives = False
        pybids_inputs = {
            "t1": {
                "wildcards": ["acquisition", "subject"],
            }
        }

        # Simplest case -- one input type, using pybids
        config = generate_inputs(
            pybids_inputs=pybids_inputs,
            bids_dir=tmpdir,
            derivatives=derivatives,
            pybids_config=str(Path(__file__).parent / "data" / "custom_config.json"),
        )
        # Order of the subjects is not deterministic
        assert config["input_lists"] == BidsListCompare({"t1": entities})
        assert config["input_zip_lists"] == {"t1": zip_list}
        assert config["input_wildcards"] == {
            "t1": {"acq": "{acq}", "subject": "{subject}"}
        }
        assert config["subj_wildcards"] == {"subject": "{subject}"}

    def test_missing_wildcards(self, tmpdir: Path):
        entities, zip_list = self.get_entities(tmpdir)

        # create config
        derivatives = False
        pybids_inputs = {
            "t1": {
                "filters": {"acquisition": "foo", "subject": "001"},
            }
        }

        # Simplest case -- one input type, using pybids
        config = generate_inputs(
            pybids_inputs=pybids_inputs,
            bids_dir=tmpdir,
            derivatives=derivatives,
            pybids_config=str(Path(__file__).parent / "data" / "custom_config.json"),
        )
        result_d = {"acq": ["foo"], "subject": ["001"]}
        assert config["input_lists"] == {"t1": {}}
        assert config["input_zip_lists"] == {"t1": {}}
        assert config["input_wildcards"] == {"t1": {}}
        assert config["subj_wildcards"] == {"subject": "{subject}"}


def test_custom_pybids_config(tmpdir: Path):
    # Generate directory
    for i in range(2):
        path = Path(
            bids(
                tmpdir, datatype="anat", subject="001", foo=str(i), suffix="T1w.nii.gz"
            )
        )
        path.parent.mkdir(parents=True, exist_ok=True)
        path.touch()

    # create config
    derivatives = False
    pybids_inputs = {
        "t1": {
            "filters": {"suffix": "T1w"},
            "wildcards": ["acquisition", "subject", "foo"],
        }
    }

    # Simplest case -- one input type, using pybids
    config = generate_inputs(
        pybids_inputs=pybids_inputs,
        bids_dir=tmpdir,
        derivatives=derivatives,
        pybids_config=(Path(__file__).parent / "data" / "custom_config.json"),
    )
    # Order of the subjects is not deterministic
    assert config["input_lists"] in [
        {"t1": {"foo": ["0", "1"], "subject": ["001"]}},
        {"t1": {"foo": ["1", "0"], "subject": ["001"]}},
    ]
    assert config["input_zip_lists"] == {
        "t1": {"foo": ["0", "1"], "subject": ["001", "001"]}
    }
    assert config["input_wildcards"] == {"t1": {"foo": "{foo}", "subject": "{subject}"}}
    # Order of the subjects is not deterministic
    assert config["subj_wildcards"] == {"subject": "{subject}"}


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
    assert config["input_lists"] == BidsListCompare(
        {"t1": {"acq": ["mprage"], "subject": ["002", "001"]}}
    )
    assert config["input_zip_lists"] == {
        "t1": {"acq": ["mprage", "mprage"], "subject": ["001", "002"]}
    }
    assert config["input_wildcards"] == {"t1": {"acq": "{acq}", "subject": "{subject}"}}
    # Order of the subjects is not deterministic
    assert set(config["subjects"]) == {"002", "001"}
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
        assert config["input_lists"] == BidsListCompare(
            {
                "t1": {"acq": ["mprage"], "subject": ["001", "002"]},
                "t2": {"subject": ["002"]},
            }
        )
        assert config["input_zip_lists"]["t1"] in [
            {"acq": ["mprage", "mprage"], "subject": ["001", "002"]},
            {"acq": ["mprage", "mprage"], "subject": ["002", "001"]},
        ]
        assert config["input_zip_lists"]["t2"] == {"subject": ["002"]}
        assert config["input_wildcards"] == {
            "t1": {"acq": "{acq}", "subject": "{subject}"},
            "t2": {"subject": "{subject}"},
        }
        # Order of the subjects is not deterministic
        assert set(config["subjects"]) == {"001", "002"}

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

        config = _get_lists_from_bids(layout, pybids_inputs)
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
        assert config["input_lists"] == BidsListCompare(
            {
                "t1": {"acq": ["mprage"], "subject": ["001", "002"]},
                "t2": {"subject": ["002"]},
            }
        )
        assert config["input_wildcards"] == {
            "t1": {"acq": "{acq}", "subject": "{subject}"},
            "t2": {"subject": "{subject}"},
        }


class TestDB:
    @pytest.fixture(autouse=True)
    def start(self, tmpdir):
        self.tmpdir = tmpdir.strpath

        # Copy over test data
        shutil.copytree("snakebids/tests/data/bids_t1w", f"{self.tmpdir}/data")
        assert filecmp.dircmp("snakebids/tests/data/bids_t1w", f"{self.tmpdir}/data")

        # Create config
        self.bids_dir = f"{self.tmpdir}/data"
        self.pybids_db = {"database_dir": "", "reset_database": False}

    def test_database_dir_blank(self):
        # Test non-saving (check db does not exist)
        _gen_bids_layout(
            bids_dir=self.bids_dir,
            derivatives=False,
            pybids_database_dir=self.pybids_db.get("database_dir"),
            pybids_reset_database=self.pybids_db.get("reset_database"),
        )
        assert not os.path.exists(self.pybids_db.get("database_dir"))

    def test_database_dir_relative(self):
        # Update config
        self.pybids_db["database_dir"] = "./.db"

        # Check to make sure db exists (relative path)
        _gen_bids_layout(
            bids_dir=self.bids_dir,
            derivatives=False,
            pybids_database_dir=self.pybids_db.get("database_dir"),
            pybids_reset_database=self.pybids_db.get("reset_database"),
        )
        assert not os.path.exists(f"{self.tmpdir}/data/.db/")

    def test_database_dir_absolute(self):
        # Update config
        self.pybids_db["database_dir"] = f"{self.tmpdir}/data/.db/"
        self.pybids_db["reset_database"] = False

        # Check to make sure db exists (absolute path)
        _gen_bids_layout(
            bids_dir=self.bids_dir,
            derivatives=False,
            pybids_database_dir=self.pybids_db.get("database_dir"),
            pybids_reset_database=self.pybids_db.get("reset_database"),
        )
        assert os.path.exists(f"{self.tmpdir}/data/.db/")

        # Test reading of old layout when changes occur
        os.makedirs(f"{self.tmpdir}/data/sub-003/anat")
        shutil.copy(
            f"{self.bids_dir}/sub-001/anat/sub-001_acq-mprage_T1w.nii.gz",
            f"{self.bids_dir}/sub-003/anat/sub-003_acq-mprage_T1w.nii.gz",
        )
        # Check to make sure new subject not cached in layout
        layout = _gen_bids_layout(
            bids_dir=self.bids_dir,
            derivatives=False,
            pybids_database_dir=self.pybids_db.get("database_dir"),
            pybids_reset_database=self.pybids_db.get("reset_database"),
        )
        assert not layout.get(subject="003")

        # Test updating of layout
        self.pybids_db["reset_database"] = True
        # Check to see if new subject in updated layout
        layout = _gen_bids_layout(
            bids_dir=self.bids_dir,
            derivatives=False,
            pybids_database_dir=self.pybids_db.get("database_dir"),
            pybids_reset_database=self.pybids_db.get("reset_database"),
        )
        assert layout.get(subject="003")
