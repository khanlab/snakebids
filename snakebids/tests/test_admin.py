import re
import sys
from argparse import ArgumentParser, Namespace
from pathlib import Path

import pytest
from hypothesis import given
from hypothesis import strategies as st
from pathvalidate import Platform, is_valid_filename
from pytest_mock.plugin import MockerFixture

from snakebids.admin import gen_parser
from snakebids.tests.helpers import allow_function_scoped


@pytest.fixture
def parser():
    return gen_parser()


class TestAdminCli:
    def test_fails_if_no_subcommand(
        self, parser: ArgumentParser, mocker: MockerFixture
    ):
        mocker.patch.object(sys, "argv", ["snakebids"])
        with pytest.raises(SystemExit):
            parser.parse_args()

    def test_fails_if_invalid_subcommand(
        self, parser: ArgumentParser, mocker: MockerFixture
    ):
        mocker.patch.object(sys, "argv", ["snakebids", "dummy"])
        with pytest.raises(SystemExit):
            parser.parse_args()

    @given(
        name=st.text()
        .filter(lambda s: not re.match(r"^[a-zA-Z_][a-zA-Z_0-9]*$", s))
        .filter(lambda s: is_valid_filename(s, Platform.LINUX))
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

    @given(name=st.text().filter(lambda s: is_valid_filename(s, Platform.LINUX)))
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

    def test_create_succeeds(self, parser: ArgumentParser, mocker: MockerFixture):
        mocker.patch.object(sys, "argv", ["snakebids", "create"])
        assert isinstance(parser.parse_args(), Namespace)

    def test_boutiques_succeeds(self, parser: ArgumentParser, mocker: MockerFixture):
        mocker.patch.object(sys, "argv", ["snakebids", "boutiques", "test.json"])
        assert isinstance(parser.parse_args(), Namespace)
