from __future__ import annotations

import copy
import logging
import subprocess
from pathlib import Path

import pytest
from pytest_mock.plugin import MockerFixture

from snakebids.app import SnakeBidsApp
from snakebids.plugins.validate import BidsValidator, InvalidBidsError
from snakebids.tests.mock.config import config


class TestBidsValidator:
    @pytest.fixture
    def app(self):
        return SnakeBidsApp(
            Path("app"),
            skip_parse_args=False,
            snakefile_path=Path("Snakefile"),
            configfile_path=Path("mock/config.yaml"),
            config=copy.deepcopy(config),
        )

    _logger = logging.getLogger(__name__)

    def test_skip_validation(self, app: SnakeBidsApp):
        # Test if validation is skipped
        app.config["plugins.validator.skip"] = True

        validator = BidsValidator()
        validator(app)
        assert "plugins.validator.success" not in app.config

    def test_validation_successful(self, app: SnakeBidsApp, mocker: MockerFixture):
        # Test successful validation
        mocker.patch("subprocess.check_call", return_value=0)

        app.config["bids_dirs"] = "path/to/bids/dir"
        app.config["plugins.validator.skip"] = False

        validator = BidsValidator()
        validator(app)
        assert app.config["plugins.validator.success"]

    def test_missing_bids_validator(
        self, app: SnakeBidsApp, mocker: MockerFixture, caplog: pytest.LogCaptureFixture
    ):
        # Test fall back to Pybids validation
        mocker.patch("subprocess.check_call", side_effect=FileNotFoundError)

        app.config["bids_dirs"] = "path/to/bids/dir"
        app.config["plugins.validator.skip"] = False

        validator = BidsValidator()
        validator(app)

        # Check validation failure
        assert not app.config["plugins.validator.success"]

        # Check logger message (should only have 1 warning message)
        assert len(caplog.records) == 1
        assert caplog.records[0].levelname == "WARNING"
        assert "Missing bids-validator installation" in caplog.records[0].message

    def test_validation_error(self, app: SnakeBidsApp, mocker: MockerFixture):
        # Test for any other bids-validation error
        mocker.patch(
            "subprocess.check_call", side_effect=subprocess.CalledProcessError(1, "")
        )

        app.config["bids_dirs"] = "path/to/bids/dir"
        app.config["plugins.validator.skip"] = False

        # Check error raised
        validator = BidsValidator()
        with pytest.raises(InvalidBidsError):
            validator(app)

        # Check validation failure
        assert not app.config["plugins.validator.success"]
