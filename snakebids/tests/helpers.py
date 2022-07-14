"""Helper functions and classes for tests
"""

from collections import UserDict
from typing import Any, Dict, Iterable, List, Tuple, Union

from snakebids import bids
from snakebids.utils.utils import BidsEntity


def get_zip_list(
    entities: Iterable[Union[BidsEntity, str]], combinations: Iterable[Tuple[str, ...]]
):
    """Return a zip list from iterables of entities and value combinations

    Parameters
    ----------
    entities : Iterable[str]
        Iterable of entities
    combinations : Iterable[Tuple[str, ...]]
        Iterable of combinations (e.g. produced using itertools.product). Each tuple
        must be the same length as entities

    Returns
    -------
    Dict[str, List[str]]
        zip_list representation of entity-value combinations
    """
    return {
        BidsEntity(entity).wildcard: list(combs)  # type: ignore
        for entity, combs in zip(entities, zip(*combinations))
    }


def setify(dic: Dict[Any, List[Any]]):
    """Convert a dict of *->lists into a dict of *->sets

    Parameters
    ----------
    dic : Dict[Any, List[Any]]
        Dict of lists to convert

    Returns
    -------
    Dict[Any, Set[Any]]
        Dict of sets
    """
    return {key: set(val) for key, val in dic.items()}


def get_bids_path(entities: Iterable[str]):
    """Get consistently ordered bids path for a group of entities

    Parameters
    ----------
    entities : Iterable[str]
        Entities to convert into a bids path

    Returns
    -------
    str
        bids path
    """

    def get_tag(entity: BidsEntity):
        # For pybids, suffixes MUST be followed by extensions, but we don't yet support
        # seperate indexing of extensions, so add a dummy extension any time there's a
        # suffix
        extension = ".foo" if entity == "suffix" else ""
        return entity.wildcard, f"{{{entity.wildcard}}}{extension}"

    return bids(
        root=".", **dict(get_tag(BidsEntity(entity)) for entity in sorted(entities))
    )


class BidsListCompare(UserDict):
    """Dict override specifically for comparing input_lists

    When comparing, all lists are converted into sets so that order doesn't matter for
    the comparison
    """

    def __eq__(self, other: Dict):
        for name, lists in other.items():
            if name not in self:
                return False
            for key, val in lists.items():
                if (
                    not isinstance(val, list)
                    or key not in self[name]
                    or set(self[name][key]) != set(val)
                ):
                    return False
        return True
