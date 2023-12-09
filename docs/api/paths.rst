Path Building
=============

.. py:currentmodule:: snakebids

.. autofunction:: snakebids.bids

.. autofunction:: bids_factory(spec)


Specs
-----

Official specs represent the evolution of the bids spec over time. Versioned :func:`bids` functions based on these specs can be imported with::

    from snakebids import bids_<version> as bids

For example, the bids function based on the ``v0_0_0`` spec can be imported using::

    from snakebids import bids_v0_0_0 as bids

The latest function can always be imported using::

    from snakebids import bids

However, this latest version is subject to breaking changes on any snakebids release, including patches. Production code should thus always import a versioned bids function.

.. automodule:: snakebids.paths.specs
    :exclude-members: latest
    :members:

.. function:: latest

    Points to the most recent spec

Types
-----

.. automodule:: snakebids
    :exclude-members: from_bids_lists, BidsComponent, BidsPartialComponent, BidsComponentRow, BidsDataset
    :members: BidsFunction

.. py:currentmodule:: snakebids.paths.utils

.. autoclass:: BidsPathEntitySpec
    :members:

.. automodule:: snakebids.paths.utils
    :noindex:
    :members: BidsPathSpec

.. class:: BidsPathSpec
