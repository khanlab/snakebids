from __future__ import annotations

from collections.abc import Iterator, Mapping
from typing import Iterable, List, Literal, TypeVar, Union, overload

import more_itertools as itx

from snakebids.types import ZipList, ZipListLike
from snakebids.utils.containers import ContainerBag, MultiSelectDict, RegexContainer

T_co = TypeVar("T_co", bound=Union[List[str], str], covariant=True)


@overload
def filter_list(
    zip_list: ZipListLike,
    filters: Mapping[str, Iterable[str] | str],
    return_indices_only: Literal[False] = ...,
    regex_search: bool = ...,
) -> ZipList:
    ...


@overload
def filter_list(
    zip_list: ZipListLike,
    filters: Mapping[str, Iterable[str] | str],
    return_indices_only: Literal[True],
    regex_search: bool = ...,
) -> list[int]:
    ...


def filter_list(
    zip_list: ZipListLike,
    filters: Mapping[str, Iterable[str] | str],
    return_indices_only: bool = False,
    regex_search: bool = False,
) -> ZipList | list[int]:
    """Filter zip_list, including only entries with provided entity values.

    Parameters
    ----------
    zip_list
        generated zip lists dict from config file to filter

    filters
        wildcard values to filter the zip lists

    return_indices_only
        return the indices of the matching wildcards

    regex_search
        Use regex matching to filter instead of the default equality check.


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
    # Save filters into memory as sets for quick access later
    if regex_search:
        filter_sets = {
            key: ContainerBag(*(RegexContainer(r) for r in itx.always_iterable(vals)))
            for key, vals in filters.items()
        }
    else:
        filter_sets = {
            key: set(itx.always_iterable(vals)) for key, vals in filters.items()
        }

    # Get a set {0,1,2,3...n-1} where n is the length of any one of the lists in
    # zip_list
    keep_indices = set(_get_zip_list_indices(zip_list)).intersection(
        *(
            {i for i, v in enumerate(zip_list[key]) if v in container}
            for key, container in filter_sets.items()
            if key in zip_list
        )
    )

    # Now we have the indices, so filter the lists
    if return_indices_only:
        return list(keep_indices)
    return MultiSelectDict(
        {key: [val[i] for i in keep_indices] for key, val in zip_list.items()}
    )


def get_filtered_ziplist_index(
    zip_list: ZipList,
    wildcards: dict[str, str],
    subj_wildcards: dict[str, str],
) -> int | list[int]:
    """Return the indices of all entries matching the filter query.

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
    subj_dict = {key: wildcards[key] for key in subj_wildcards}

    # now filter the list based on subj_wildcards
    zip_list_filtered = filter_list(zip_list, subj_dict)

    # get the index of the wildcard from this filtered list
    indices = filter_list(zip_list_filtered, wildcards, return_indices_only=True)
    if len(indices) == 1:
        return indices[0]
    return indices


def _get_zip_list_indices(zip_list: ZipListLike) -> Iterator[int]:
    """Convert a zip_list into its indices.

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
    zip_list : dict[str, list[str]]
        Zip_list to be converted

    Yields
    ------
    integers
        The indices of the zip_list
    """
    if not zip_list:
        return

    # Retrieve the first zip_list (all zip_lists should be the same, so it doesn't
    # matter which one we take)
    sample_zip_list = zip_list[itx.first(zip_list)]

    # Return a sequence of numbers 0, 1, 2, 3, ... n-1 where n is the length of the
    # zip list
    yield from range(len(sample_zip_list))
