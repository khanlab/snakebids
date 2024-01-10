BIDS App Bootstrapping
======================


.. automodule:: snakebids.app
    :members:

BIDS App
--------

.. currentmodule:: snakebids.bidsapp

.. module:: snakebids.bidsapp

.. autofunction:: app

.. autoclass:: snakebids.bidsapp.run._Runner()

    .. current module:: snakebids.bidsapp

    .. autoattribute:: pm
    .. autoattribute:: parser
    .. autoattribute:: config
    .. autoattribute:: argument_groups

    Action Methods
    --------------

    Plugins are only run when calling the action methods. These methods trigger
    running of the plugins up to a specified point. For example,
    :meth:`~_Runner.build_parser` runs the :func:`~hookspecs.initialize_config`
    and :func:`~hookspecs.add_cli_arguments` hooks. It can thus be used to build
    the parser without actually parsing any arguments.

    Plugins are always run in the same order. Action methods will always trigger
    the entire plugin chain up until their stopping point. Importantly, plugins
    will only ever be run once, even if action methods are called multiple
    times. For example, if :meth:`~_Runner.build_parser` is called, then
    :meth:`~_Runner.parse_args`, the parser will only be built once.

    .. automethod:: build_parser
    .. automethod:: parse_args
    .. automethod:: run

