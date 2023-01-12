# flake8: noqa
__submodules__ = ["construct_bids", "filtering", "input_generation", "datasets"]

__ignore__ = ["T_co"]

# <AUTOGEN_INIT>
from snakebids.core.construct_bids import bids, print_boilerplate
from snakebids.core.datasets import BidsComponent, BidsDataset, BidsDatasetDict
from snakebids.core.filtering import filter_list, get_filtered_ziplist_index
from snakebids.core.input_generation import (
    generate_inputs,
    get_wildcard_constraints,
    write_derivative_json,
)

__all__ = [
    "BidsComponent",
    "BidsDataset",
    "BidsDatasetDict",
    "bids",
    "filter_list",
    "generate_inputs",
    "get_filtered_ziplist_index",
    "get_wildcard_constraints",
    "print_boilerplate",
    "write_derivative_json",
]

# </AUTOGEN_INIT>
