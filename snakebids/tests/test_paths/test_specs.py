from snakebids.paths import specs
from snakebids.paths.config import set_bids_spec
from snakebids.paths.presets import bids
from snakebids.paths.specs import v0_0_0
from snakebids.paths.utils import find_entity


def test_all_entries_define_entity():
    spec = v0_0_0()
    for item in spec:
        assert "entity" in item


def test_subject_dir_can_be_excluded():
    spec = v0_0_0(subject_dir=False)
    subject = find_entity(spec, "subject")
    assert subject.get("dir") is False


def test_session_dir_can_be_excluded():
    spec = v0_0_0(session_dir=False)
    session = find_entity(spec, "session")
    assert session.get("dir") is False


def test_spec_can_be_set_with_str():
    set_bids_spec("v0_0_0")
    assert bids(acquisition="foo") == "acquisition-foo"
    set_bids_spec("v0_10_1")
    assert bids(acquisition="foo") == "acq-foo"


def test_spec_can_be_set_with_obj():
    set_bids_spec(specs.v0_0_0())
    assert bids(acquisition="foo") == "acquisition-foo"
    set_bids_spec(specs.v0_10_1())
    assert bids(acquisition="foo") == "acq-foo"
