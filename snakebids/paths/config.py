from __future__ import annotations

from typing import Literal, TypedDict, cast

from typing_extensions import TypeAlias

from snakebids.paths import specs
from snakebids.paths.factory import BidsFunction, bids_factory
from snakebids.paths.utils import BidsPathSpec

# <AUTOUPDATE>
# The code between these tags is automatically generated. Do not
# manually edit
# To update, run::
#
#       poetry run poe update-bids
#

VALID_SPECS: TypeAlias = Literal["v0_10_1", "v0_0_0", "latest"]


# </AUTOUPDATE>
class _Config(TypedDict):
    active_spec: BidsPathSpec
    bids_func: BidsFunction


_latest_spec = specs.latest()
_config = _Config(active_spec=_latest_spec, bids_func=bids_factory(_latest_spec))


def set_bids_spec(spec: BidsPathSpec | VALID_SPECS):
    """Set the spec to be used by path generation functions (such as :func:`~snakebids.bids`).

    Parameters
    ----------
    spec
        Either a spec object, or the name of a builtin :ref:`spec <specs>`
    """  # noqa: E501
    if isinstance(spec, str):
        spec = cast("BidsPathSpec", getattr(specs, spec)())
    _config["active_spec"] = spec
    _config["bids_func"] = bids_factory(spec)


def get_bids_spec() -> BidsPathSpec:
    """Get the currently active BIDS spec."""
    return _config["active_spec"]


def get_bids_func() -> BidsFunction:
    """Get the currently active BIDS path build function."""
    return _config["bids_func"]
