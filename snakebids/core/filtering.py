import itertools as it
import operator as op
import re
from typing import Dict, List, Tuple, Union, overload

import more_itertools as itx
from typing_extensions import Literal

from snakebids.utils.utils import matches_any


@overload
def filter_list(
    zip_list,
    filters: Union[Dict[str, str], Dict[str, List[str]]],
    return_indices_only: Literal[False] = ...,
    regex_search: bool = ...,
) -> Dict[str, Tuple[str]]:
    ...


@overload
def filter_list(
    zip_list,
    filters: Union[Dict[str, str], Dict[str, List[str]]],
    return_indices_only: Literal[True] = ...,
    regex_search: bool = ...,
) -> Tuple[int]:
    ...


def filter_list(
    zip_list: Dict[str, List[str]],
    filters: Union[Dict[str, str], Dict[str, List[str]]],
    return_indices_only=False,
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
        ...     'dir': ('AP', 'PA', 'AP', 'PA'),
        ...     'acq': ('98', '98', '99', '99'),
        ...     'subject': ('01', '01', '01', '01')
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
        ...     'dir': ('AP', 'PA', 'AP', 'PA'),
        ...     'acq': ('98', '98', '98', '98'),
        ...     'subject': ('01', '01', '02', '02')
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
        ...     'dir': ('AP', 'AP', 'AP', 'AP'),
        ...     'acq': ('98', '98', '99', '99'),
        ...     'subject': ('01', '02', '01', '02')
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
        ...     'dir': (),
        ...     'acq': (),
        ...     'subject': ()
        ... }
        True
    """

    # Unzip values to group them by path instead of by wildcard type:
    #   ['AP', 'PA', 'AP', 'PA']            ['AP', '01']
    #   ['01', '01', '02', 'PA']    ->      ['PA', '01']
    #                                       ['AP', '02']
    #                                       ['PA', '02']
    value_groups = zip(*zip_list.values())

    # Select a match function based on the value of regex_search
    if regex_search:
        match_func = re.match
    else:
        match_func = op.eq

    def filter_values(args: Tuple[int, Tuple[str]]):
        # This function is used to filter each value group obtained above. The first
        # member of args is an enumeration integer not used for filtering.
        value_group = args[1]

        # Reassociate each member of the value group with its wildcard key
        for key, value in zip(zip_list.keys(), value_group):
            # If the key has a filter, check if the value matches any of the filters
            if key in filters and not matches_any(
                value, [*itx.always_iterable(filters[key])], match_func
            ):
                # If not, immediately return False
                return False
        # If every value in the group passes the filters, return True
        return True

    # Run the above filter function. Enumerate is used to track which indices were kept
    filtered = [*filter(filter_values, enumerate(value_groups))]
    if len(filtered):
        # Unzip filtered to seperated the indices from the filtered_values
        filtered_values: Tuple[Tuple[str]]
        indices, filtered_values = zip(*filtered)  # type: ignore
        if return_indices_only:
            return indices
        # Unzip filtered_values to restore to the original structure. These are zipped
        # with keys and turned back into a dict
        return dict(zip(zip_list.keys(), zip(*filtered_values)))

    # If no values survivied filtering, we return a dict with all the keys mapping to
    # empty tuples
    return dict(it.zip_longest(zip_list.keys(), [], fillvalue=tuple()))


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
