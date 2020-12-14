Installation
============

Template snakebids BIDS App 


Requirements
------------

Docker (Mac/Windows/Linux) or Singularity (Linux)

Docker:
^^^^^^^

Pull the container:

.. code-block::

   docker pull khanlab/app_name:latest

do a dry run, printing the command at each step:

.. code-block::

   docker run -it --rm -v PATH_TO_BIDS_DIR:/bids:ro -v PATH_TO_OUTPUT_DIR:/output khanlab/app_name:latest /bids /output participant -np 

run it with maximum number of cores:

.. code-block::

   docker run -it --rm -v PATH_TO_BIDS_DIR:/bids:ro -v PATH_TO_OUTPUT_DIR:/output khanlab/app_name:latest /bids /output participant -p --cores all


Singularity:
^^^^^^^^^^^^

Pull the container:

.. code-block::

   singularity pull khanlab_app_name_latest.sif docker://khanlab/app_name:latest

do a dry run, printing the command at each step:

.. code-block::

   singularity run -e khanlab_app_name_latest.sif khanlab/app_name:latest PATH_TO_BIDS_DIR PATH_TO_OUTPUT_DIR participant -np 

run it with maximum number of cores:

.. code-block::

   singularity run -e khanlab_app_name_latest.sif khanlab/app_name:latest PATH_TO_BIDS_DIR PATH_TO_OUTPUT_DIR participant  -p --cores all



