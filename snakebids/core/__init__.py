__submodules__ = ["filtering", "input_generation", "datasets"]

__ignore__ = ["T_co"]

from snakebids.core.datasets import (
    BidsComponent,
    BidsComponentRow,
    BidsDataset,
    BidsDatasetDict,
    BidsPartialComponent,
)

# <AUTOGEN_INIT>
from snakebids.core.filtering import (
    filter_list,
    get_filtered_ziplist_index,
)
from snakebids.core.input_generation import (
    generate_inputs,
    get_wildcard_constraints,
    write_derivative_json,
)

__all__ = [
    "BidsComponent",
    "BidsComponentRow",
    "BidsDataset",
    "BidsDatasetDict",
    "BidsPartialComponent",
    "filter_list",
    "generate_inputs",
    "get_filtered_ziplist_index",
    "get_wildcard_constraints",
    "write_derivative_json",
]
# </AUTOGEN_INIT>
