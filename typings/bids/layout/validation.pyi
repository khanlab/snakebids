"""
This type stub file was generated by pyright.
"""

"""Functionality related to validation of BIDSLayouts and BIDS projects."""
MANDATORY_BIDS_FIELDS = ...
MANDATORY_DERIVATIVES_FIELDS = ...
EXAMPLE_BIDS_DESCRIPTION = ...
EXAMPLE_DERIVATIVES_DESCRIPTION = ...
DEFAULT_LOCATIONS_TO_IGNORE = ...

def absolute_path_deprecation_warning(): ...
def indexer_arg_deprecation_warning(): ...
def validate_root(root, validate): ...
def validate_derivative_paths(paths, layout=..., **kwargs): ...
def validate_indexing_args(ignore, force_index, root): ...
