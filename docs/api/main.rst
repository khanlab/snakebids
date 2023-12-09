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
    app


Path Building
-------------

.. Need to manually create this table because bids does not have a proper docstring

===================================  ================================
:func:`bids <bids>`                  Generate bids or bids-like paths
:func:`bids_factory <bids_factory>`  Create new :func:`bids` functions according to a spec
===================================  ================================


Dataset Creation
----------------

.. autosummary::
    :nosignatures:

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



BIDS App Booststrapping
-----------------------

.. currentmodule:: snakebids.app

.. autosummary::
    :nosignatures:
    :recursive:

    SnakeBidsApp
