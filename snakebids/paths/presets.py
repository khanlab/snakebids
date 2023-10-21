from __future__ import annotations

import functools as ft
from typing import TYPE_CHECKING

from snakebids.paths import specs
from snakebids.paths.factory import bids_factory
from snakebids.paths.specs import latest

# <AUTOUPDATE>
# The code between these tags is automatically generated. Do not
# manually edit
# To update, run::
#
#       poetry run poe update_bids
#

if not TYPE_CHECKING:
    __all__ = [  # noqa:F822
        "bids_v0_0_0",
        "bids",
    ]

    def __dir__():
        return __all__


# </AUTOUPDATE>


@ft.lru_cache
def __getattr__(name: str):
    if name == "bids":
        return bids_factory(latest())
    prefix = name[:5]
    version = name[5:]
    if prefix != "bids_" or (spec := getattr(specs, version, None)) is None:
        err = f"module '{__name__}' has no attribute '{name}'"
        raise AttributeError(err)
    return bids_factory(spec())
