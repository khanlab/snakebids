Data Structures
===============

.. py:currentmodule:: snakebids

.. autoclass:: BidsComponent
    :members:
    :exclude-members: input_wildcards, input_lists, input_name, input_path, input_zip_lists
    :inherited-members:

.. dropdown:: Legacy ``BidsComponents`` properties
    :icon: info
    :class-title: sd-outline-info

    The following properties are historical aliases of ``BidsComponents`` properties. There are no current plans to deprecate them, but new code should avoid them.

    .. autoproperty:: BidsComponent.input_zip_lists

    .. autoproperty:: BidsComponent.input_wildcards

    .. autoproperty:: BidsComponent.input_name

    .. autoproperty:: BidsComponent.input_path

    .. autoproperty:: BidsComponent.input_lists


.. autoclass:: BidsPartialComponent

.. autoclass:: BidsComponentRow
    :members:
    :exclude-members: zip_lists

.. autoclass:: BidsDataset
    :members:
    :exclude-members: input_wildcards, input_lists, input_path, input_zip_lists

.. dropdown:: Legacy ``BidsDataset`` properties
    :icon: info
    :class-title: sd-outline-info

    The following properties are historical aliases of :class:`~snakebids.BidsDataset` properties. There are no current plans to deprecate them, but new code should avoid them.

    .. autoproperty:: BidsDataset.input_zip_lists

    .. autoproperty:: BidsDataset.input_wildcards

    .. autoproperty:: BidsDataset.input_path

    .. autoproperty:: BidsDataset.input_lists

.. autoclass:: BidsDatasetDict
    :members:
