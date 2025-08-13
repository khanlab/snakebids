from __future__ import annotations

import argparse
import re
import sys
from argparse import ArgumentParser, Namespace
from pathlib import Path

import more_itertools as itx
import pytest
from hypothesis import given
from hypothesis import strategies as st
from pathvalidate import Platform, is_valid_filename
from pytest_mock.plugin import MockerFixture

from snakebids.admin import create_descriptor, gen_parser
from snakebids.tests.helpers import allow_function_scoped


@pytest.fixture
def parser():
    return gen_parser()


def test_fails_if_no_subcommand(parser: ArgumentParser, mocker: MockerFixture):
    mocker.patch.object(sys, "argv", ["snakebids"])
    with pytest.raises(SystemExit):
        parser.parse_args()


def test_fails_if_invalid_subcommand(parser: ArgumentParser, mocker: MockerFixture):
    mocker.patch.object(sys, "argv", ["snakebids", "dummy"])
    with pytest.raises(SystemExit):
        parser.parse_args()


class TestCreateCommand:
    @given(
        name=st.text()
        .filter(lambda s: not re.match(r"^[a-zA-Z_][a-zA-Z_0-9]*$", s))
        .filter(lambda s: is_valid_filename(s, Platform.LINUX))
        .filter(lambda s: s not in {".", ".."})
    )
    @allow_function_scoped
    def test_create_fails_with_invalid_filename(
        self,
        parser: ArgumentParser,
        mocker: MockerFixture,
        name: str,
        tmp_path: Path,
        capsys: pytest.CaptureFixture[str],
    ):
        mocker.patch.object(sys, "argv", ["snakebids", "create", str(tmp_path / name)])
        args = parser.parse_args()
        with pytest.raises(SystemExit):
            args.func(args)
        capture = capsys.readouterr()
        assert "valid python module" in capture.err
        assert name in capture.err

    @given(
        name=st.text()
        .filter(lambda s: is_valid_filename(s, Platform.LINUX))
        .filter(lambda s: s not in {".", ".."})
    )
    @allow_function_scoped
    def test_create_fails_missing_parent_dir(
        self,
        parser: ArgumentParser,
        mocker: MockerFixture,
        name: str,
        tmp_path: Path,
        capsys: pytest.CaptureFixture[str],
    ):
        path = tmp_path / name / "sub"
        mocker.patch.object(sys, "argv", ["snakebids", "create", str(path)])
        args = parser.parse_args()
        with pytest.raises(SystemExit):
            args.func(args)
        capture = capsys.readouterr()
        assert "does not exist" in capture.err
        assert str(path.parent) in capture.err

    def test_create_fails_when_markers_in_snakebids_version(
        self,
        parser: ArgumentParser,
        mocker: MockerFixture,
        capsys: pytest.CaptureFixture[str],
    ):
        mocker.patch.object(
            sys,
            "argv",
            ["snakebids", "create", "name", "--snakebids-version", "... ; markers"],
        )
        args = parser.parse_args()
        with pytest.raises(SystemExit):
            args.func(args)
        capture = capsys.readouterr()
        assert "may not specify markers" in capture.err

    def test_create_fails_when_get_spec_has_at_sign(
        self,
        parser: ArgumentParser,
        mocker: MockerFixture,
        capsys: pytest.CaptureFixture[str],
    ):
        mocker.patch.object(
            sys,
            "argv",
            ["snakebids", "create", "name", "--snakebids-version", "@ git+...@...@..."],
        )
        args = parser.parse_args()
        with pytest.raises(SystemExit):
            args.func(args)
        capture = capsys.readouterr()
        assert "in git requirement" in capture.err

    def test_create_fails_when_snakebids_version_specifies_extras(
        self,
        parser: ArgumentParser,
        mocker: MockerFixture,
        capsys: pytest.CaptureFixture[str],
    ):
        mocker.patch.object(
            sys,
            "argv",
            ["snakebids", "create", "name", "--snakebids-version", "[...] ..."],
        )
        args = parser.parse_args()
        with pytest.raises(SystemExit):
            args.func(args)
        capture = capsys.readouterr()
        assert "may not specify extras" in capture.err

    @given(
        name=st.from_regex(r"^[a-zA-Z_][a-zA-Z_0-9]*$"),
        version=st.text(st.characters(blacklist_characters=["@", ";", "["])).filter(
            lambda s: not s.startswith("-")
        )
        | st.none(),
    )
    @allow_function_scoped
    def test_create_calls_copier_correctly(
        self,
        parser: ArgumentParser,
        mocker: MockerFixture,
        name: str,
        version: str | None,
        tmp_path: Path,
        capsys: pytest.CaptureFixture[str],
    ):
        import snakebids
        from snakebids.admin import copier

        run_copy = mocker.patch.object(copier, "run_copy")
        path = tmp_path / name
        args = ["snakebids", "create", str(path)]
        if version is not None:
            args.extend(["--snakebids-version", version])
        mocker.patch.object(sys, "argv", args)
        args = parser.parse_args()
        args.func(args)
        capture = capsys.readouterr()
        assert "Creating Snakebids app at" in capture.err
        run_copy.assert_called_once_with(
            str(Path(itx.first(snakebids.__path__), "project_template")),
            path,
            data={
                "app_full_name": name,
                **({"snakebids_version": version} if version is not None else {}),
            },
            unsafe=True,
        )

    def test_create_succeeds(self, parser: ArgumentParser, mocker: MockerFixture):
        mocker.patch.object(sys, "argv", ["snakebids", "create"])
        assert isinstance(parser.parse_args(), Namespace)

    def test_boutiques_succeeds(self, parser: ArgumentParser, mocker: MockerFixture):
        mocker.patch.object(sys, "argv", ["snakebids", "boutiques", "test.json"])
        assert isinstance(parser.parse_args(), Namespace)


def test_boutiques_descriptor():
    with pytest.raises(SystemExit):
        create_descriptor(argparse.Namespace())
