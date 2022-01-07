# flake8: noqa
__submodules__ = ["bids", "filtering", "input_generation"]

# <AUTOGEN_INIT>
from snakebids.core.bids import bids, print_boilerplate
from snakebids.core.filtering import filter_list, get_filtered_ziplist_index
from snakebids.core.input_generation import (
    generate_inputs,
    get_wildcard_constraints,
    write_derivative_json,
)

__all__ = [
    "bids",
    "filter_list",
    "generate_inputs",
    "get_filtered_ziplist_index",
    "get_wildcard_constraints",
    "print_boilerplate",
    "write_derivative_json",
]

# </AUTOGEN_INIT>
