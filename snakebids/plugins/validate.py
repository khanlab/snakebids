import json
import logging
import subprocess
import tempfile

from snakebids.app import SnakeBidsApp

_logger = logging.getLogger(__name__)


class InvalidBidsError(Exception):
    """Error raised if an input BIDS dataset is invalid."""


class BidsValidator:
    """Perform BIDS validation of dataset, falling back to the pybids
    version of validation if the node-version of bids-validator is
    not found.

    Parameters
    ----------
    app
        Snakebids application to be run
    """

    def __init__(self) -> None:
        pass

    def __call__(self, app: SnakeBidsApp) -> None:
        self.app = app

        # Skip bids validation
        if self.app.config["skip_bids_validation"]:
            return

        try:
            validator_config_dict = {"ignoredFiles": ["/participants.tsv"]}

            with tempfile.NamedTemporaryFile(mode="w+", suffix=".json") as temp:
                temp.write(json.dumps(validator_config_dict))
                temp.flush()

            subprocess.check_call(
                ["bids-validator", self.app.config["bids_dirs"], "-c", temp.name]
            )

            # If successfully bids-validation performed, skip pybids validation
            self.app.config["bids_validator_success"] = True

        # If the bids-validator call can't be made
        except FileNotFoundError:
            self.app.config["bids_validator_success"] = False
            _logger.warning(
                "Missing bids-validator installation - falling back to pybids "
                "validation."
            )
        # Any other bids-validator error
        except subprocess.CalledProcessError as err:
            self.app.config["bids_validator_success"] = False
            raise InvalidBidsError from err
