import warnings

import pytest
from pytest_mock import MockerFixture

from snakebids.paths import bids_factory, specs
from snakebids.paths._config import reset_bids_spec, set_bids_spec
from snakebids.paths._presets import bids
from snakebids.paths._utils import find_entity
from snakebids.paths.specs import v0_0_0


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


def test_using_include_subject_dir_raises_warning():
    with pytest.warns(UserWarning, match="include_session_dir and include_subject_dir"):
        bids(subject="001", include_subject_dir=False)
    with pytest.warns(UserWarning, match="include_session_dir and include_subject_dir"):
        bids(session="001", include_session_dir=False)


class TestCustomEntityWarnings:
    def test_using_custom_entities_with_default_bids_raises_warning(self):
        reset_bids_spec()
        with pytest.warns(UserWarning, match="spec has not been explicitly declared"):
            bids(foo="bar")

    def test_no_warning_when_spec_declared(self):
        reset_bids_spec()
        set_bids_spec("v0_0_0")
        with warnings.catch_warnings():
            warnings.simplefilter("error")
            bids(foo="bar")

    def test_no_warning_when_bids_explicitly_generated(self):
        reset_bids_spec()
        with warnings.catch_warnings():
            warnings.simplefilter("error")
            bids_factory(specs.v0_0_0())(foo="bar")

    def test_no_warning_in_interactive_mode(self, mocker: MockerFixture):
        reset_bids_spec()
        mocker.patch(
            "snakebids.paths._factory.in_interactive_session", return_value=True
        )
        with warnings.catch_warnings():
            warnings.simplefilter("error")
            bids(foo="bar")
