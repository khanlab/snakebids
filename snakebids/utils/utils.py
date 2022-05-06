import functools as ft
import importlib.resources
import json
from typing import Any, Callable, Dict, Iterable


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


def matches_any(
    item: Any,
    match_list: Iterable[Any],
    match_func: Callable[[Any, Any], Any],
    *args: Any
):
    for match in match_list:
        if match_func(match, item, *args):
            return True
    return False


def get_match_search_func(
    match_list: Iterable[Any], match_func: Callable[[Any, Any], Any]
):
    """Return a match search function suitable for use in filter

    Parameters
    ----------
    match_list : list
        list of items to search for matches
    match_func : callable
        Any callable that takes two args and returns truthy/falsy values. The first arg
        will be drawn from match_list, the second will be the value being matched

    Returns
    -------
    callable
        Takes as a single arg a value to be matched. It will be compared to every item
        in match list using match_func
    """
    match_list = list(match_list)

    def inner(item: Any):
        return matches_any(item, match_list, match_func)

    return inner
