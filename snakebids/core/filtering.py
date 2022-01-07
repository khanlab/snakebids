def filter_list(zip_list, wildcards, return_indices_only=False):
    """This function is used when you are expanding over some subset of the
    wildcards i.e. if your output file doesn't contain all the wildcards in
    input_wildcards

    Parameters
    ----------
    zip_list : dict
        generated zip lists dict from config file to filter

    wildcards : dict
        wildcard values to filter the zip lists

    return_indices_only : bool, default=False
        return the indices of the matching wildcards

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
    keep_indices = []
    for key, val in wildcards.items():
        # get indices where wildcard exists
        if key not in zip_list.keys():
            continue
        indices = [i for i, v in enumerate(zip_list[key]) if v == val]
        if len(keep_indices) == 0:
            keep_indices = indices
        else:
            keep_indices = [x for x in indices if x in keep_indices]
    # Now we have the indices, so filter the lists
    if return_indices_only:
        return keep_indices
    return {
        key: [zip_list[key][i] for i in keep_indices] for key, val in zip_list.items()
    }


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
