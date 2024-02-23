================
Internals
================

.. note::

    These types are mostly used internally. The API most users will need is documented in :ref:`the main API page<main api>`, but these are listed here for reference (and some of the main API items link here).


utils
-----

.. automodule:: snakebids.utils.utils
    :members:


snakemake_io
------------

.. automodule:: snakebids.utils.snakemake_io
    :members:


exceptions
----------

.. automodule:: snakebids.exceptions
    :members:

types
-----

.. automodule:: snakebids.types
    :exclude-members: ZipList, InputsConfig, ZipListLike
    :members:


.. Workaround to get links to the TypeAliases working.
    As documented here: https://github.com/sphinx-doc/sphinx/issues/10785

    The issue is that autodoc does not automatically create the links to the
    aliases correctly, because they get documented as :py:data:. So we document
    the three of them as classes to give autodoc something to link to, then use
    css to hide those entries. And the actual visible documentation is generated
    below with the :noindex: tag to prevent overlap with the pseudo-class
    documentation

.. automodule:: snakebids.types
    :members: InputsConfig, ZipList, ZipListLike
    :noindex:

.. class:: ZipList

.. class:: InputsConfig

.. class:: ZipListLike
