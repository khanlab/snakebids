"""Helper functions and classes for tests
"""

from collections import UserDict
from typing import Dict


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
