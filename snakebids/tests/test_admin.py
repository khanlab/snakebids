from __future__ import absolute_import

import sys
from argparse import ArgumentParser

import pytest
from configargparse import Namespace
from pytest_mock.plugin import MockerFixture

from snakebids.admin import gen_parser


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

    def test_create_succeeds(self, parser: ArgumentParser, mocker: MockerFixture):
        mocker.patch.object(sys, "argv", ["snakebids", "create"])
        assert isinstance(parser.parse_args(), Namespace)

    def test_boutiques_succeeds(self, parser: ArgumentParser, mocker: MockerFixture):
        mocker.patch.object(sys, "argv", ["snakebids", "boutiques", "test.json"])
        assert isinstance(parser.parse_args(), Namespace)
