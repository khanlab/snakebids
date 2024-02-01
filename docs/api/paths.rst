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

BIDS specs control the formatting of paths produced by the :func:`~snakebids.bids` function. They specify the order of recognized entities, placing ``ses-X`` after ``sub-Y``, for instance, no matter what order they are specified in the function. Unrecognized entitites are placed in the order specified in the function call.

Because of this, each addition of entities to the spec presents a potentially breaking change. Suppose an entity called ``foo`` were added to the spec. Calls to :func:`~snakebids.bids` with ``foo`` as an argument would place the entity at the end of the path:

.. code-block:: python

    from snakebids import bids

    # Before foo is in the spec
    bids(
        subject="001",
        session="1",
        label="WM",
        foo="bar",
        suffix="data.nii.gz",
    ) == "sub-001_ses-1_label-WM_foo-bar_data.nii.gz"

The addition of ``foo`` to the spec might move the position of the entity forward in the output:

.. code-block:: python

    # After foo is in the spec
    bids(
        subject="001",
        session="1",
        label="WM",
        foo="bar",
        suffix="data.nii.gz",
    ) == "sub-001_ses-1_foo-bar_label-WM_data.nii.gz"

This would change the output paths of workflow using this function call, causing a breaking change in workflow behaviour.

To ensure stable path generation across releases, Snakebids ships with versioned specs that can be explicitly set using :func:`snakebids.set_bids_specs`. These specs are named after the snakebids version they release with. By default, :func:`~snakebids.bids` will always use the latest spec, but production code should generally declare the spec to be used by the workflow:

.. code-block:: python

    from snakebids import set_bids_spec
    set_bids_spec("v0_0_0")

This is especially true of workflows using custom entities. To emphasize this, a warning is issued in python scripts and apps using such entities without declaring a spec version.


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
