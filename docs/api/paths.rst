Path Building
=============

.. py:currentmodule:: snakebids

.. autofunction:: snakebids.bids

.. autofunction:: bids_factory(spec)

.. autofunction:: set_bids_spec

.. _specs:

Specs
-----

.. py:currentmodule:: snakebids.paths.specs

Official specs represent the evolution of the BIDS spec over time. Currently, :func:`~snakebids.bids` uses the :func:`v0_0_0` by default, however, in a future release, this will change to always use the latest spec. Specs may be updated at any time, without warning, even on patch releases, so it's important for production code to explicitly specify the BIDS spec version in use. This is done using :func:`~snakebids.set_bids_spec`:

.. code-block:: python

    from snakebids import set_bids_spec
    set_bids_spec("v0_0_0")


.. automodule:: snakebids.paths.specs
    :exclude-members: latest
    :members:

.. function:: latest

    Points to the most recent spec

Types
-----

.. py:currentmodule:: snakebids

.. autoclass:: BidsFunction
    :members:

.. py:currentmodule:: snakebids.paths

.. autoclass:: BidsPathEntitySpec
    :members:

.. automodule:: snakebids.paths
    :noindex:
    :members: BidsPathSpec

.. class:: BidsPathSpec
