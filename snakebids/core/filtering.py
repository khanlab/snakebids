import operator as op
import re
from typing import Dict, List, TypeVar, Union, overload

import more_itertools as itx
from typing_extensions import Literal

from snakebids.utils.utils import matches_any

# pylint: disable=invalid-name
T_co = TypeVar("T_co", bound=Union[List[str], str], covariant=True)


@overload
def filter_list(
    zip_list,
    filters: Dict[str, T_co],
    return_indices_only: Literal[False] = ...,
    regex_search: bool = ...,
) -> Dict[str, List[str]]:
    ...


@overload
def filter_list(
    zip_list,
    filters: Dict[str, T_co],
    return_indices_only: Literal[True] = ...,
    regex_search: bool = ...,
) -> List[int]:
    ...


def filter_list(
    zip_list: Dict[str, List[str]],
    filters: Dict[str, T_co],
    return_indices_only: bool = False,
    regex_search=False,
):
    """This function is used when you are expanding over some subset of the
    wildcards i.e. if your output file doesn't contain all the wildcards in
    input_wildcards

    Parameters
    ----------
    zip_list : dict
        generated zip lists dict from config file to filter

    filters : dict
        wildcard values to filter the zip lists

    return_indices_only : bool, default=False
        return the indices of the matching wildcards

    regex_search : bool, default=False
        Use regex matching to filter instead of the default equality check.

    Returns
    -------
    dict
        zip list with non-matching elements removed

    Examples
    --------
    ::

        >>> import snakebids

    Filtering to get all ``subject='01'`` scans::

        >>> snakebids.filter_list(
        ...     {
        ...         'dir': ['AP','PA','AP','PA', 'AP','PA','AP','PA'],
        ...         'acq': ['98','98','98','98','99','99','99','99'],
        ...         'subject': ['01','01','02','02','01','01','02','02' ]
        ...     },
        ...     {'subject': '01'}
        ... ) == {
        ...     'dir': ['AP', 'PA', 'AP', 'PA'],
        ...     'acq': ['98', '98', '99', '99'],
        ...     'subject': ['01', '01', '01', '01']
        ... }
        True

    Filtering to get all ``acq='98'`` scans::

        >>> snakebids.filter_list(
        ...     {
        ...         'dir': ['AP','PA','AP','PA', 'AP','PA','AP','PA'],
        ...         'acq': ['98','98','98','98','99','99','99','99'],
        ...         'subject': ['01','01','02','02','01','01','02','02' ]
        ...     },
        ...     {'acq': '98'}
        ... ) == {
        ...     'dir': ['AP', 'PA', 'AP', 'PA'],
        ...     'acq': ['98', '98', '98', '98'],
        ...     'subject': ['01', '01', '02', '02']
        ... }
        True

    Filtering to get all ``dir=='AP'`` scans::

        >>> snakebids.filter_list(
        ...     {
        ...         'dir': ['AP','PA','AP','PA', 'AP','PA','AP','PA'],
        ...         'acq': ['98','98','98','98','99','99','99','99'],
        ...         'subject': ['01','01','02','02','01','01','02','02' ]
        ...     },
        ...     {'dir': 'AP'}
        ... ) == {
        ...     'dir': ['AP', 'AP', 'AP', 'AP'],
        ...     'acq': ['98', '98', '99', '99'],
        ...     'subject': ['01', '02', '01', '02']
        ... }
        True

    Filtering to get all ``subject='03'`` scans (i.e. no matches)::

        >>> snakebids.filter_list(
        ...     {
        ...         'dir': ['AP','PA','AP','PA', 'AP','PA','AP','PA'],
        ...         'acq': ['98','98','98','98','99','99','99','99'],
        ...         'subject': ['01','01','02','02','01','01','02','02' ]
        ...     },
        ...     {'subject': '03'}
        ... ) == {
        ...     'dir': [],
        ...     'acq': [],
        ...     'subject': []
        ... }
        True
    """
    if regex_search:
        match_func = re.match
    else:
        match_func = op.eq

    keep_indices = set.intersection(
        # Get a set {0,1,2,3...n-1} where n is the length of any one of the lists in
        # zip_list
        {*_get_zip_list_indices(zip_list)},
        *(
            {
                i
                for i, v in enumerate(zip_list[key])
                if matches_any(v, itx.always_iterable(val), match_func)
            }
            for key, val in filters.items()
            if key in zip_list
        )
    )

    # Now we have the indices, so filter the lists
    if return_indices_only:
        return list(keep_indices)
    return {key: [val[i] for i in keep_indices] for key, val in zip_list.items()}


def get_filtered_ziplist_index(zip_list, wildcards, subj_wildcards):
    """Use this function when you have wildcards for a single scan instance,
    and want to know the index of that scan, amongst that subject's scan
    instances.

    Parameters
    ----------
    zip_list : dict
        lists for scans in a dataset, zipped to get each instance
    wildcards: dict
        wildcards for the single instance for querying it's index
    subj_wildcards: dict
        keys of this dictionary are used to pick out the subject/(session)
        from the wildcards

    Examples
    --------

    >>> import snakebids

    In this example, we have a dataset where with scans from two subjects,
    where each subject has ``dir-AP`` and ``dir-PA`` scans, along with
    ``acq-98`` and ``acq-99``:

    * ``sub-01_acq-98_dir-AP_dwi.nii.gz``
    * ``sub-01_acq-98_dir-PA_dwi.nii.gz``
    * ``sub-01_acq-99_dir-AP_dwi.nii.gz``
    * ``sub-01_acq-99_dir-PA_dwi.nii.gz``
    * ``sub-02_acq-98_dir-AP_dwi.nii.gz``
    * ``sub-02_acq-98_dir-PA_dwi.nii.gz``
    * ``sub-02_acq-99_dir-AP_dwi.nii.gz``
    * ``sub-02_acq-99_dir-PA_dwi.nii.gz``

    The ``zip_list`` produced by ``generate_inputs()`` is the set of entities
    that when zipped together, e.g. with ``expand(path, zip, **zip_list)``,
    produces the entity combinations that refer to each scan::

        {
            'dir': ['AP','PA','AP','PA', 'AP','PA','AP','PA'],
            'acq': ['98','98','98','98','99','99','99','99'],
            'subject': ['01','01','02','02','01','01','02','02']
        }

    The ``filter_list()`` function produces a subset of the entity
    combinations as a filtered zip list. This is used e.g. to get all the
    scans for a single subject.

    This ``get_filtered_ziplist_index()`` function performs ``filter_list()``
    twice:

    1. Using the ``subj_wildcards`` (e.g.: ``'subject': '{subject}'``) to get
       a subject/session-specific ``zip_list``.
    2. To return the indices from that list of the matching wildcards.

    In this example, if the wildcards parameter was::

        {'dir': 'PA', 'acq': '99', 'subject': '01'}

    Then the first (subject/session-specific) filtered list provides this zip
    list::

        {
            'dir': ['AP','PA','AP','PA'],
            'acq': ['98','98','99','99'],
            'subject': ['01','01','01','01']
        }

    which has 4 combinations, and thus are indexed from 0 to 3.

    The returned value would then be the index (or indices) that matches the
    wildcards. In this case, since the wildcards were
    ``{'dir': 'PA', 'acq': '99', 'subject':'01'}``, the return index is 3. ::

        >>> snakebids.get_filtered_ziplist_index(
        ...     {
        ...         'dir': ['AP','PA','AP','PA', 'AP','PA','AP','PA'],
        ...         'acq': ['98','98','98','98','99','99','99','99'],
        ...         'subject': ['01','01','02','02','01','01','02','02' ]
        ...     },
        ...     {'dir': 'PA', 'acq': '99', 'subject': '01'},
        ...     {'subject': '{subject}' }
        ... )
        3
    """
    # get the subject/(session) dict:
    subj_dict = {key: wildcards[key] for key in subj_wildcards.keys()}

    # now filter the list based on subj_wildcards
    zip_list_filtered = filter_list(zip_list, subj_dict)

    # get the index of the wildcard from this filtered list
    indices = filter_list(zip_list_filtered, wildcards, return_indices_only=True)
    if len(indices) == 1:
        return indices[0]
    return indices


def _get_zip_list_indices(zip_list: Dict[str, List[str]]):
    """Convert a zip_list into its indices

    Generates a sequence of numbers from 0 up to the length of the zip_lists. For
    example, the zip list:

    zip_list = {
        "subject": ["001", "002", "001", "002"],
        "contrast": ["T1w", "T1w", "T2w", "T2w"]
    }

    would lead to:

    (0, 1, 2, 3)

    Parameters
    ----------
    zip_list : Dict[str, List[str]]
        Zip_list to be converted

    Yields
    -------
    integers
        The indices of the zip_list
    """

    if not zip_list:
        return

    # Retrieve the first zip_list (all zip_lists should be the same, so it doesn't
    # matter which one we take)
    # pylint: disable=stop-iteration-return
    sample_zip_list = zip_list[next(iter(zip_list))]

    # Return a sequence of numbers 0, 1, 2, 3, ... n-1 where n is the length of the
    # zip list
    yield from range(len(sample_zip_list))
