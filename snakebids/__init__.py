# type: ignore
__version__ = "0.0.0"

__submodules__ = ["core", "paths"]

# isort: split
# <AUTOGEN_INIT>
import lazy_loader

__getattr__, __dir__, __all__ = lazy_loader.attach_stub(__name__, __file__)

__all__ = [
    "BidsComponent",
    "BidsComponentRow",
    "BidsDataset",
    "BidsDatasetDict",
    "BidsPartialComponent",
    "bids",
    "bids_v0_0_0",
    "filter_list",
    "generate_inputs",
    "get_filtered_ziplist_index",
    "get_wildcard_constraints",
    "write_derivative_json",
]
# </AUTOGEN_INIT>
