Overview
========

The ``bids`` function generates a BIDS-like filepath corresponding to its keyword arguments. The generated filepath has the form::

    [root]/[sub-{subject}]/[ses-{session]/[prefix]_[sub-{subject}]_[ses-{session}]_[{key}-{val}_ ... ]_[suffix]

Use cases of the ``bids`` function include, at simplest, replacing a hard-coded BIDS file with an invocation of the ``bids`` function, so ::

    "data/sub-01/ses-01/func/sub-01_ses-01_task-rest_acq-01_run-1_bold.nii.gz"

could become ::

    bids(root="data", subject="01", "session"="01", datatype="func", task="rest", acq="01", run="1", suffix="bold.nii.gz")

If you wanted to specify that a rule should run on a BIDS file from any subject, session, acquisition, task, and run, you could change those keyword arguments to be snakemake wildcards::

    bids(root="data", subject="{subject}", session="{session}", datatype="func", task="{task}", acq="{acq}", run="{run}", suffix="bold.nii.gz")

Using the subject and session keywords as wildcards is common enough that snakebids pre-populates a config variable (``subj_wildcards``) with these wildcards, allowing the ``bids`` call to look like the following::

    bids(root="data", datatype="func", task="{task}", acq="{acq}", run="{run}", suffix="bold.nii.gz", **config["subj_wildcards"])

Now if you want to process all inputs of a given form regardless of how their wildcards resolve, snakebids can provide the necessary arguments the base Snakemake ``expand`` function based on what's grabbed from the input dataset. The required information is in the ``input_zip_lists`` config variable. As an example, to specify the output of a rule that preprocesses BOLD images (as specified in the example configuration), the following would resolve the subject, session, acquisition, task, and run wildcards by expanding over ``input_zip_lists``::

    expand(bids(root="output", datatype="func", desc="preproc", acq="{acq}", task="{task}", run="{run}", **config["subj_wildcards"]), zip, config["input_zip_lists"]["bold"])
