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
    "BidsFunction",
    "BidsPartialComponent",
    "bids",
    "bids_factory",
    "filter_list",
    "generate_inputs",
    "get_bids_func",
    "get_bids_spec",
    "get_filtered_ziplist_index",
    "get_wildcard_constraints",
    "set_bids_spec",
    "write_derivative_json",
]
# </AUTOGEN_INIT>
