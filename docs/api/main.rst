.. _main api:

================
API
================


snakebids
---------

.. autoclass:: snakebids.BidsComponent
    :members:
    :exclude-members: input_wildcards, input_lists, input_name, input_path, input_zip_lists
    :inherited-members:

.. dropdown:: Legacy ``BidsComponents`` properties
    :icon: info
    :class-title: sd-outline-info

    The following properties are historical aliases of ``BidsComponents`` properties. There are no current plans to deprecate them, but new code should avoid them.

    .. autoproperty:: snakebids.BidsComponent.input_zip_lists

    .. autoproperty:: snakebids.BidsComponent.input_wildcards

    .. autoproperty:: snakebids.BidsComponent.input_name

    .. autoproperty:: snakebids.BidsComponent.input_path

    .. autoproperty:: snakebids.BidsComponent.input_lists


.. autoclass:: snakebids.BidsPartialComponent

.. autoclass:: snakebids.BidsComponentRow
    :members:
    :exclude-members: zip_lists

.. autoclass:: snakebids.BidsDataset
    :members:
    :exclude-members: input_wildcards, input_lists, input_path, input_zip_lists

.. dropdown:: Legacy ``BidsDataset`` properties
    :icon: info
    :class-title: sd-outline-info

    The following properties are historical aliases of :class:`~snakebids.BidsDataset` properties. There are no current plans to deprecate them, but new code should avoid them.

    .. autoproperty:: snakebids.BidsDataset.input_zip_lists

    .. autoproperty:: snakebids.BidsDataset.input_wildcards

    .. autoproperty:: snakebids.BidsDataset.input_path

    .. autoproperty:: snakebids.BidsDataset.input_lists


.. automodule:: snakebids
    :exclude-members: from_bids_lists, BidsComponent, BidsPartialComponent, BidsComponentRow, BidsDataset
    :members:


app
---

.. automodule:: snakebids.app
    :members:
