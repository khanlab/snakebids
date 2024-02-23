Plugins
=======

.. note::

    This page consists of plugins distributed with Snakebids. Externally developed and distributed plugins may be missing from this page.

Core
----
These plugins provide essential bids app features and will be needed by most apps.

.. currentmodule:: snakebids.plugins

.. autoclass:: BidsArgs
    :members:

.. autoclass:: ComponentEdit
    :members:

.. autoclass:: CliConfig
    :members:

.. autoclass:: Pybidsdb
    :members:


Utility
-------
These plugins add additional features that may be helpful to most bids apps.

.. autoclass:: BidsValidator
    :members:


Workflow Integrations
---------------------
These plugins enable integrations with specific workflow managers.

.. autoclass:: SnakemakeBidsApp
    :members:


Plugin Development
------------------
.. function:: snakebids.bidsapp.hookimpl()

    Marker to be imported and used in plugins (and for own implementations).

    See :func:`pluggy.HookimplMarker` for more information. Its parameters are
    represented here fore reference.

    :param optionalhook:
        If ``True``, a missing matching hook specification will not result
        in an error (by default it is an error if no matching spec is
        found). See :ref:`optionalhook`.

    :param tryfirst:
        If ``True``, this hook implementation will run as early as possible
        in the chain of N hook implementations for a specification. See
        :ref:`callorder`.

    :param trylast:
        If ``True``, this hook implementation will run as late as possible
        in the chain of N hook implementations for a specification. See
        :ref:`callorder`.

    :param wrapper:
        If ``True`` ("new-style hook wrapper"), the hook implementation
        needs to execute exactly one ``yield``. The code before the
        ``yield`` is run early before any non-hook-wrapper function is run.
        The code after the ``yield`` is run after all non-hook-wrapper
        functions have run. The ``yield`` receives the result value of the
        inner calls, or raises the exception of inner calls (including
        earlier hook wrapper calls). The return value of the function
        becomes the return value of the hook, and a raised exception becomes
        the exception of the hook. See :external:ref:`hookwrappers`.

    :param hookwrapper:
        If ``True`` ("old-style hook wrapper"), the hook implementation
        needs to execute exactly one ``yield``. The code before the
        ``yield`` is run early before any non-hook-wrapper function is run.
        The code after the ``yield`` is run after all non-hook-wrapper
        function have run  The ``yield`` receives a :class:`pluggy.Result` object
        representing the exception or result outcome of the inner calls
        (including earlier hook wrapper calls). This option is mutually
        exclusive with ``wrapper``. See :external:ref:`old_style_hookwrappers`.

    :param specname:
        If provided, the given name will be used instead of the function
        name when matching this hook implementation to a hook specification
        during registration. See :ref:`specname`.


Specs
~~~~~
The spec contains the hooks recognized by snakebids. Each spec has a specific function name and a set of available arguments, and will be called at a specific time during app initialization.

.. automodule:: snakebids.bidsapp.hookspecs
    :members:
