import json
import logging
import subprocess
import tempfile

from snakebids.app import SnakeBidsApp

_logger = logging.getLogger(__name__)


class InvalidBidsError(Exception):
    """Error raised if an input BIDS dataset is invalid."""


def bids_validate(app: SnakeBidsApp, bids_dir: str) -> None:
    """Perform validation of dataset. Initial attempt at validation performed
    with node-version of bids-validator. If not found, will fall back to Python
    version of validation (same as pybids).

    Parameters
    ----------
    app
        Snakebids application to be run
    bids_dir
        BIDS organized directory to be validated
    """

    # Skip bids validation
    if app.config["skip_bids_validation"]:
        return

    try:
        validator_config_dict = {"ignoredFiles": ["/participants.tsv"]}

        with tempfile.NamedTemporaryFile(mode="w+", suffix=".json") as temp:
            temp.write(json.dumps(validator_config_dict))
            temp.flush()

        subprocess.check_call(["bids-validator", str(bids_dir), "-c", temp.name])
    # If the bids-validator call can't be made
    except FileNotFoundError:
        _logger.warning(
            "Bids-validator does not appear to be installed - will use python "
            "validation."
        )
    # Any other bids-validator error
    except subprocess.CalledProcessError as err:
        raise InvalidBidsError from err
