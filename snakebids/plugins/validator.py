from __future__ import annotations

import json
import logging
import subprocess as sp
import tempfile
from typing import Any

import attr

from snakebids import bidsapp
from snakebids.exceptions import SnakebidsPluginError

_logger = logging.getLogger(__name__)


class InvalidBidsError(SnakebidsPluginError):
    """Error raised if an input BIDS dataset is invalid."""


@attr.define
class BidsValidator:
    """Perform validation of a BIDS dataset using the bids-validator.

    If the dataset is not valid according to the BIDS specifications, an
    InvalidBidsError is raised.

    Parameters
    ----------
    raise_invalid_bids : bool
        Flag to indicate whether InvalidBidsError should be raised if BIDS
        validation fails. Default to True.

    """

    raise_invalid_bids: bool = attr.field(default=True)

    def __eq__(self, other: Any):
        return isinstance(other, self.__class__)

    @bidsapp.hookimpl
    def finalize_config(self, config: dict[str, Any]) -> None:
        """Perform BIDS validation of dataset.

        Raises
        ------
        InvalidBidsError
            Raised when the input BIDS directory does not pass validation with
            the bids-validator
        """
        # Skip bids validation
        if config["plugins.validator.skip"]:
            return

        validator_config_dict = {"ignoredFiles": ["/participants.tsv"]}

        with tempfile.NamedTemporaryFile(mode="w+", suffix=".json") as temp:
            temp.write(json.dumps(validator_config_dict))
            temp.flush()
            try:
                sp.check_call(["bids-validator", config["bids_dir"], "-c", temp.name])

                # If successfully bids-validation performed
                config["plugins.validator.success"] = True
            except FileNotFoundError:
                # If the bids-validator call can't be made
                config["plugins.validator.success"] = False
                _logger.warning(
                    "Missing bids-validator installation - falling back to pybids "
                    "validation."
                )
            # Any other bids-validator error
            except sp.CalledProcessError as err:
                config["plugins.validator.success"] = False
                if self.raise_invalid_bids:
                    raise InvalidBidsError from err
