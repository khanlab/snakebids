from __future__ import annotations

import argparse
import subprocess as sp
from pathlib import Path

import pytest
from pytest_mock.plugin import MockerFixture

from snakebids.plugins.validator import BidsValidator, InvalidBidsError


def test_skip_validation():
    val = BidsValidator()
    parser = argparse.ArgumentParser()
    val.add_cli_arguments(parser)
    parsed = parser.parse_args(["--skip-bids-validation"])
    assert parsed.__dict__["plugins.validator.skip"] is True


class TestFinalizeConfig:
    def test_skip_validation(self):
        val = BidsValidator()
        # Test if validation is skipped
        config = {"plugins.validator.skip": True}

        val.finalize_config(config)
        assert "plugins.validator.success" not in config

    def test_validation_successful(self, mocker: MockerFixture):
        # Test successful validation
        mocker.patch("subprocess.check_call", return_value=0)
        val = BidsValidator()
        # Test if validation is skipped
        config = {"plugins.validator.skip": False, "bids_dir": Path()}

        val.finalize_config(config)

        assert config["plugins.validator.success"]

    def test_missing_bids_validator(
        self, mocker: MockerFixture, caplog: pytest.LogCaptureFixture
    ):
        # Test fall back to Pybids validation
        mocker.patch("subprocess.check_call", side_effect=FileNotFoundError)
        val = BidsValidator()
        # Test if validation is skipped
        config = {"plugins.validator.skip": False, "bids_dir": Path()}

        val.finalize_config(config)

        assert not config["plugins.validator.success"]

        # Check logger message (should only have 1 warning message)
        assert len(caplog.records) == 1
        assert caplog.records[0].levelname == "WARNING"
        assert "Missing bids-validator installation" in caplog.records[0].message

    def test_raise_validation_error(self, mocker: MockerFixture):
        # Test for any other bids-validation error
        mocker.patch("subprocess.check_call", side_effect=sp.CalledProcessError(1, ""))
        val = BidsValidator()
        # Test if validation is skipped
        config = {"plugins.validator.skip": False, "bids_dir": Path()}

        with pytest.raises(InvalidBidsError):
            val.finalize_config(config)

        # Check validation failure
        assert not config["plugins.validator.success"]

    def test_ignore_validation_error(self, mocker: MockerFixture):
        # Test for any other bids-validation error
        mocker.patch("subprocess.check_call", side_effect=sp.CalledProcessError(1, ""))
        val = BidsValidator(raise_invalid_bids=False)
        # Test if validation is skipped
        config = {"plugins.validator.skip": False, "bids_dir": Path()}
        val.finalize_config(config)

        # Check validation failure
        assert not config["plugins.validator.success"]
