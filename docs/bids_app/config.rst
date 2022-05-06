Configuration
=============

Snakebids is configured with a YAML (or JSON) file that extends the standard `snakemake config file <https://snakemake.readthedocs.io/en/stable/snakefiles/configuration.html#standard-configuration>`_ with variables that snakebids uses to parse an input BIDS dataset and expose the snakebids workflow to the command line.

Config Variables
----------------

pybids_inputs: A dictionary that describes each type of input you want to grab from an input BIDS dataset. The value of each item should be a dictionary with keys ``filters`` and ``wildcards``. The value of ``filters`` should be a dictionary where each key corresponds to a BIDS entity, and the value specifies which values of that entity should be grabbed. The value of ``wildcards`` should be a list of BIDS entities that should be present, but specifies that all values of that entity should be grabbed. The dictionary for each input is passed directly to `PyBIDS' get() function <https://bids-standard.github.io/pybids/generated/bids.layout.BIDSLayout.html#bids.layout.BIDSLayout.get>`_.

In the following (YAML-formatted) example, the ``bold`` input type is specified. BIDS files with the datatype ``func``, suffix ``bold``, and extension ``.nii.gz`` will be grabbed, and the ``subject``, ``session``, ``acquisition``, ``task``, and ``run`` entities of those files will be left as wildcards. ::

    bold:
      filters:
        suffix: 'bold'
        extension: '.nii.gz'
        datatype: 'func'
      wildcards:
        - subject
        - session
        - acquisition
        - task
        - run

pybids_db_dir: PyBIDS allows for the use of a cached layout to be used in order to reduce the time required to index a BIDS dataset. A path (if provided) to save the ``PyBIDS`` layout. If ``None`` or ``''`` is provided, the layout is not saved / used. The path provided must be absolute, otherwise the database will not be used. Note, this is a variable used for an opt-in feature and must first be uncommented in the ``snakebids.yml`` file.

pybids_db_reset: A boolean determining whether the existing layout should be be updated. Default behaviour does not update the existing database if one is used. Note, this is a variable used for an opt-in feature and must first be uncommented in the ``snakebids.yml`` file.

analysis_levels: A list of analysis levels in the BIDS app. Typically, this will include participant and/or group. Note that the default (YAML) configuration file expects this mapping to be identified with the anchor ``analysis_levels`` to be aliased by ``parse_args``.

targets_by_analysis_level: A mapping from the name of each ``analysis_level`` to the list of rules or files to be run for that analysis level.

parse_args: A dictionary of command-line parameters to make available as part of the BIDS app. Each item of the mapping is passed to `argparse's add_argument function <https://docs.python.org/3/library/argparse.html#the-add-argument-method>`_. A number of default entries are present in a new snakebids project's config file that structure the BIDS app's CLI, but additional command-line arguments can be added as necessary.

debug: A boolean that determines whether debug statements are printed during parsing. Should be disabled (False) if you're generating DAG visualization with snakemake.

derivatives: A boolean (or path(s) to derivatives datasets) that determines whether snakebids will search in the derivatives subdirectory of the input dataset.
