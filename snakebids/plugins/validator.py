import json
import logging
import subprocess
import tempfile

from snakebids.app import SnakeBidsApp
from snakebids.exceptions import SnakebidsPluginError

_logger = logging.getLogger(__name__)


class InvalidBidsError(SnakebidsPluginError):
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

    def __init__(self, raise_invalid_bids: bool = True) -> None:
        self.raise_invalid_bids = raise_invalid_bids

    def __call__(self, app: SnakeBidsApp) -> None:
        # Skip bids validation
        if app.config["plugins.validator.skip"]:
            return

        validator_config_dict = {"ignoredFiles": ["/participants.tsv"]}

        with tempfile.NamedTemporaryFile(mode="w+", suffix=".json") as temp:
            temp.write(json.dumps(validator_config_dict))
            temp.flush()
            try:
                subprocess.check_call(
                    ["bids-validator", app.config["bids_dirs"], "-c", temp.name]
                )

                # If successfully bids-validation performed
                app.config["plugins.validator.success"] = True
            except FileNotFoundError:
                # If the bids-validator call can't be made
                app.config["plugins.validator.success"] = False
                _logger.warning(
                    "Missing bids-validator installation - falling back to pybids "
                    "validation."
                )
            # Any other bids-validator error
            except subprocess.CalledProcessError as err:
                app.config["plugins.validator.success"] = False
                if self.raise_invalid_bids:
                    raise InvalidBidsError from err
