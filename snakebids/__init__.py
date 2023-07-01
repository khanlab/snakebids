__version__ = "0.0.0"

__submodules__ = ["core", "paths"]
__explicit__ = ["factory_specs"]

from snakebids.paths import specs as factory_specs

# isort: split
# <AUTOGEN_INIT>
from snakebids.core import (
    BidsComponent,
    BidsComponentRow,
    BidsDataset,
    BidsDatasetDict,
    BidsPartialComponent,
    filter_list,
    generate_inputs,
    get_filtered_ziplist_index,
    get_wildcard_constraints,
    write_derivative_json,
)
from snakebids.paths import bids

__all__ = [
    "BidsComponent",
    "BidsComponentRow",
    "BidsDataset",
    "BidsDatasetDict",
    "BidsPartialComponent",
    "bids",
    "factory_specs",
    "filter_list",
    "generate_inputs",
    "get_filtered_ziplist_index",
    "get_wildcard_constraints",
    "write_derivative_json",
]

# </AUTOGEN_INIT>
