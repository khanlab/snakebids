from __future__ import annotations

import functools as ft
from typing import TYPE_CHECKING

from snakebids.paths._templates import spec_func
from snakebids.paths._utils import BidsPathSpec, find_entity, get_spec_path, load_spec

# <AUTOUPDATE>
# The code between these tags is automatically generated. Do not
# manually edit
# To update, run::
#
#       poetry run poe update-bids
#
if not TYPE_CHECKING:
    __all__ = [  # noqa:F822
        "v0_0_0",
        "v0_11_0",
        "latest",
        "LATEST",
    ]

    def __dir__():
        return __all__


_SPECS = ["v0_0_0", "v0_11_0"]
# LATEST = "v0_11_0"
# </AUTOUPDATE>

# To automatically use latest spec as "LATEST", remove this line and uncomment the
# generated line in scripts/update_bids.py
LATEST = "v0_0_0"


@ft.lru_cache
def __getattr__(name: str):
    """Allow dynamic retrieval of latest spec."""
    if name == "latest":
        name = LATEST

    if name not in _SPECS:
        msg = f"module '{__name__}' has no attribute '{name}'"
        raise AttributeError(msg)

    spec_config = load_spec(get_spec_path(name))

    spec = spec_config["spec"]

    def get_spec(subject_dir: bool = True, session_dir: bool = True) -> BidsPathSpec:
        if not subject_dir:
            find_entity(spec, "subject")["dir"] = False

        if not session_dir:
            find_entity(spec, "session")["dir"] = False

        return spec

    get_spec.__doc__ = spec_func.format_doc(spec_config)
    get_spec.__name__ = name

    return get_spec
