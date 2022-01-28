import functools as ft
import importlib.resources
import json
from typing import Dict


@ft.lru_cache(None)
def read_bids_tags(bids_json=None) -> Dict[str, str]:
    """Read the bids tags we are aware of from a JSON file. This is used
    specifically for compatibility with pybids, since some tag keys are
    different from how they appear in the file name, e.g. ``subject`` for
    ``sub``, and ``acquisition`` for ``acq``.

    Parameters
    ----------
    bids_json : str, optional
        Path to JSON file to use, if not specified will use
        ``bids_tags.json`` in the snakebids module.

    Returns
    -------
    dict:
        Dictionary of bids tags"""
    if bids_json:
        with bids_json.open("r") as infile:
            bids_tags = json.load(infile)
        return bids_tags
    with importlib.resources.open_text("snakebids.resources", "bids_tags.json") as file:
        bids_tags = json.load(file)
    return bids_tags
