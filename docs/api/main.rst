.. _main api:


================
API
================

.. py:currentmodule:: snakebids

.. toctree::
    :hidden:

    paths
    creation
    manipulation
    structures


Path Creation
-------------

.. Need to manually create this table because bids does not have a proper docstring

.. ===================================  ================================
.. :func:`bids <bids>`                  Generate bids or bids-like paths
.. :func:`bids_factory <bids_factory>`  Create new :func:`bids` functions according to a spec
.. ===================================  ================================

.. autosummary::
    :nosignatures:

    bids
    bids_factory

Dataset Creation
----------------

.. autosummary::

    generate_inputs

Dataset Manipulation
--------------------

.. autosummary::
    :nosignatures:

    filter_list
    get_filtered_ziplist_index

Data Structures
---------------

.. autosummary::
    :nosignatures:

    BidsComponent
    BidsPartialComponent
    BidsComponentRow
    BidsDataset
    BidsDatasetDict



app
---

.. automodule:: snakebids.app
    :members:
