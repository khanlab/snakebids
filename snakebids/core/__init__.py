# type: ignore
__submodules__ = ["filtering", "input_generation", "datasets"]

__ignore__ = ["T_co"]


# <AUTOGEN_INIT>
import lazy_loader

__getattr__, __dir__, __all__ = lazy_loader.attach_stub(__name__, __file__)

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
