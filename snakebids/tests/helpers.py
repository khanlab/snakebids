"""Helper functions and classes for tests
"""

from collections import UserDict
from typing import Any, Dict, Iterable, List, Tuple


def get_zip_list(entities: Iterable[str], combinations: Iterable[Tuple[str, ...]]):
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
    return {entity: list(combs) for entity, combs in zip(entities, zip(*combinations))}


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
