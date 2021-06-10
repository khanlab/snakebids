Tutorial
========

In this example we will make a workflow to smooth `bold` scans from a bids dataset.

We will start by creating a simple rule, then make this more generalizable in each step. To begin with, this is the command we are using to smooth a bold scan. ::

    fslmaths ../bids/sub-001/func/sub-001_task-rest_run-1_bold.nii.gz -s 2.12 results/sub-001/func/sub-001_task-rest_run-1_fwhm-5mm_bold.nii.gz



Step 0:
-------
In this rule, we start by creating a rule that is effectively hard-coding the paths for input and output to re-create the command as above.

.. literalinclude:: step0/Snakefile
  :language: python

.. asciinema:: step0/step0.cast

