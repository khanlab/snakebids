from __future__ import annotations

import warnings
from pathlib import Path

import pytest
from pytest_mock import MockerFixture

import snakebids.paths._templates.spec_func as spec_template
from snakebids.paths import bids_factory, specs
from snakebids.paths._config import reset_bids_spec, set_bids_spec
from snakebids.paths._presets import bids
from snakebids.paths._utils import find_entity, get_spec_path, load_spec
from snakebids.paths.specs import v0_0_0


@pytest.fixture(autouse=True)
def _spec_reset():  # type: ignore
    reset_bids_spec()


def test_reset_bids_spec_clears_cache(mocker: MockerFixture):
    load_spec_spy = mocker.spy(specs, "load_spec")
    _ = specs.v0_0_0
    load_spec_spy.assert_called_once()
    load_spec_spy.reset_mock()
    _ = specs.v0_0_0
    load_spec_spy.assert_not_called()
    reset_bids_spec()
    load_spec_spy.reset_mock()
    _ = specs.v0_0_0
    load_spec_spy.assert_called_once()


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
    set_bids_spec("v0_11_0")
    assert bids(acquisition="foo") == "acq-foo"


def test_spec_can_be_set_with_obj():
    set_bids_spec(specs.v0_0_0())
    assert bids(acquisition="foo") == "acquisition-foo"
    set_bids_spec(specs.v0_11_0())
    assert bids(acquisition="foo") == "acq-foo"


def test_using_include_subject_dir_raises_warning():
    with pytest.warns(UserWarning, match="include_session_dir and include_subject_dir"):
        bids(subject="001", include_subject_dir=False)
    with pytest.warns(UserWarning, match="include_session_dir and include_subject_dir"):
        bids(session="001", include_session_dir=False)


def test_include_subject_dir_can_remove_dir():
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        assert len(Path(bids(subject="001", include_subject_dir=False)).parents) == 1
        assert len(Path(bids(session="001", include_session_dir=False)).parents) == 1


class TestSpecDocstrings:
    def clean(self, item: str):
        return " ".join(item.split())

    def specfuncs(self):
        for attr in dir(specs):
            print(attr)
            if attr.startswith("v") or attr == "latest":
                yield getattr(specs, attr)

    def test_all_specs_have_docstrings(self):
        for spec in self.specfuncs():
            assert isinstance(spec.__doc__, str)

    def test_all_specs_have_format_example(self):
        for spec in self.specfuncs():
            assert "Formatted as::" in spec.__doc__

    def test_spec_has_correct_docstring(self):
        spec = load_spec(get_spec_path("v0_0_0"))
        assert (
            self.clean(specs.v0_0_0.__doc__).startswith(self.clean(spec["description"]))  # type: ignore
        )

    def test_spec_with_description_gets_default_docstring(self, mocker: MockerFixture):
        spec = load_spec(get_spec_path("v0_0_0"))
        _ = spec.pop("description", None)
        mocker.patch.object(specs, "load_spec", return_value=spec)
        assert self.clean(specs.v0_0_0.__doc__).startswith(  # type: ignore
            self.clean(spec_template.DEFAULT_DESCRIPTION.format(version="v0.0.0"))
        )

    def test_no_format_example_when_no_docstring_parser(self, mocker: MockerFixture):
        mocker.patch.object(
            spec_template, "_import_docstring_parser", side_effect=ImportError()
        )
        # Need to reset after mocking because the reset collects the latest spec
        reset_bids_spec()
        for spec in self.specfuncs():
            assert "Formatted as::" not in spec.__doc__, spec


class TestCustomEntityWarnings:
    def test_using_custom_entities_with_default_bids_raises_warning(self):
        with pytest.warns(UserWarning, match="spec has not been explicitly declared"):
            bids(foo="bar")

    def test_no_warning_when_spec_declared(self):
        set_bids_spec("v0_0_0")
        with warnings.catch_warnings():
            warnings.simplefilter("error")
            bids(foo="bar")

    def test_no_warning_when_bids_explicitly_generated(self):
        with warnings.catch_warnings():
            warnings.simplefilter("error")
            bids_factory(specs.v0_0_0())(foo="bar")

    def test_no_warning_in_interactive_mode(self, mocker: MockerFixture):
        mocker.patch(
            "snakebids.paths._factory.in_interactive_session", return_value=True
        )
        with warnings.catch_warnings():
            warnings.simplefilter("error")
            bids(foo="bar")
