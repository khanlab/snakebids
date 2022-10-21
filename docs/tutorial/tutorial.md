Getting started
===============

In this example we will make a workflow to smooth ``bold`` scans from a bids dataset.

We will start by creating a simple rule, then make this more generalizable in each step. To begin with, this is the command we are using to smooth a bold scan.

    fslmaths ../bids/sub-001/func/sub-001_task-rest_run-1_bold.nii.gz -s 2.12 results/sub-001/func/sub-001_task-rest_run-1_fwhm-5mm_bold.nii.gz

This command performs smoothing with a sigma=2.12 Gaussian kernel (equivalent to 5mm FWHM, with sigma=fwhm/2.355), and saves the smoothed file as ``results/sub-001/func/sub-001_task-rest_run-1_fwhm-5mm_bold.nii.gz``.

Prerequisites
-------------

- To go through this tutorial you only need to have snakemake and snakebids installed (``pip install snakebids``)
- Since we will just make use of snakemake dry-runs (``-n`` option ), the tutorial data are just placeholders (zero-sized files), and you don't actually need to have FSL installed for ``fslmaths``.

*  TODO: instructions to get or recreate the bids tutorial data

Part I: Snakemake
=================

Step 0: a basic non-generic workflow
------------------------------------

In this rule, we start by creating a rule that is effectively hard-coding the paths for input and output to re-create the command as above.

In this rule we have an ``input:`` section for input **files**, a ``params:`` section for **non-file**  parameters, and an ``output:`` section for output files. The shell section is used to build the shell command, and can refer to the input, output or params using curly braces. Note that any of these can also be named inputs, but we have only used this for the ``sigma`` parameter in this case.

```{literalinclude} step0/Snakefile
  :language: python
```

With this rule in our Snakefile, we can then run ``snakemake -np`` to execute a dry-run of the workflow. Here the ``-n`` specifies dry-run, and the ``-p`` prints any shell commands that are to be executed.

```{asciinema} step0/step0.cast
```

When we invoke ``snakemake``, it uses the first rule in the snakefile as the ``target`` rule. The target rule is what snakemake uses as a starting point to determine what other rules need to be run in order to generate the inputs. We'll learn a bit more about that in the next step.

So far, we just have a fancy way of specifying the exact same command we started with, so there is no added benefit (yet). But we will soon add to this rule to make it more generalizable.

Step 1: adding wildcards
------------------------

First step to make the workflow generalizeable is to replace the hard-coded identifiers (e.g. the subject, task and run) with wildcards.

In the Snakefile, we can replace sub-001 with sub-{subject}, and so forth for task and run.  Now the rule is generic for any subject, task, or run.

```{literalinclude} step1/Snakefile
  :language: python
```

However, now, if we try to execute (dry-run) the workflow as before, we get an error. This is because the ``target`` rule now has wildcards in it. So snakemake is unable to determine what rules need to be run to generate the inputs, since the wildcards can take any value.

```{asciinema} step1/step1.cast
```

So for the time being, we will make use of the snakemake command-line argument to specify ``targets``, and specify the file we want generated from the command-line, by running::

    snakemake -np results/sub-001/func/sub-001_task-rest_run-1_fwhm-5mm_bold.nii.gz

We can now even try running this for another subject by changing the target file. ::

    snakemake -np results/sub-002/func/sub-002_task-rest_run-1_fwhm-5mm_bold.nii.gz

Try using a subject that doesn't exist in our bids dataset, what happens?

Now, try changing the output smoothing value, e.g. ``fwhm-10mm``, and see what happens.
As expected the command still uses a smoothing value of 2.12, since that has been hard-coded, but we will see how to rectify this in the next step.

Step 2: adding a params function
--------------------------------

As we noted, the sigma parameter needs to be computed from the FWHM. We can use a function to do this. Functions can be used for any ``input`` or ``params``, and must take ``wildcards`` as an input argument, which provides a mechanism to pass the wildcards (determined from the output file) to the function.

We can thus define a simple function that returns a string representing ``FWHM/2.355`` as follows:

```py
def calc_sigma_from_fwhm(wildcards):
    return f'{float(wildcards.fwhm)/2.355:0.2f}'
```

Note 1: We now have to make the fwhm in the output filename a wildcard, so that it can be passed to the function (via the wildcards object).

Note 2: We have to convert the fwhm to float, since all wildcards are always strings (since they are parsed from the output filename).

Once we have this function, we can replace the hardcoded ``2.12`` with the name of the function::

```py
    params:
        sigma = calc_sigma_from_fwhm
```

Here is the full Snakefile:

```{literalinclude} step2/Snakefile
  :language: python
```

Now try running the workflow again, with fwhm-5 as well as fwhm-10.

```{asciinema} step2/step2.cast
```

Step 3: adding a target rule
----------------------------

Now we have a generic rule, but it is pretty tedious to have to type out the filename of each target from the command-line in order to use it.

This is where target rules come in. If you recall from earlier, the first rule in a workflow is interpreted as the target rule, so we just need to add a dummy rule to the Snakefile that has all the target files as inputs. It is a dummy rule since it doesn't have any outputs or any command to run itself, but snakemake will take these input files, and determine if any other rules in the workflow can generate them (considering any wildcards too).

In this case, we have a BIDS dataset with two runs (run-1, run-2), and suppose we wanted to compute smoothing with several different FWHM kernels (5,10,15,20). We can thus make a target rule that has all these resulting filenames as inputs.

A very useful function in snakemake is ``expand()`` (`see Snakemake docs <https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#the-expand-function>`_). It is a way to perform array expansion to create lists of strings (input filenames).

    rule all:
        input:
            expand(
              'results/sub-{subject}/func/sub-{subject}_task-{task}_run-{run}_fwhm-{fwhm}mm_bold.nii.gz',
              subject='001',
              task='rest',
              run=[1,2],
              fwhm=[5,10,15,20]
            )


Now, we don't need to specify any targets from the command-line, and can just run::

    snakemake -np

```{asciinema} step3/step3.cast
```

The entire Snakefile for reference is:

```{literalinclude} step3/Snakefile
  :language: python
```

(step_4)=
Step 4: adding a config file
----------------------------

We have a functional workflow, but suppose you need to configure or run it on another bids dataset with different subjects, tasks, runs, or you want to run it for different smoothing values. You have to actually modify your workflow in order to do this.

It is a better practice instead to keep your configuration variables separate from the actual workflow. Snakemake supports this by allowing for a separate config file (can be YAML or JSON, here we will use YAML), where we can store any dataset specific configuration. Then to apply it for a new purpose, you can simply update the config file.

To do this, we simply add a line to our workflow::

    configfile: 'config.yml'

Snakemake will then handle reading it in, and making the configuration variables available via dictionary called ``config``.

In our config file, we will add variables for everything in the target rule ``expand()``::

    subjects:
      - '001'

    tasks:
      - rest

    runs:
      - 1
      - 2

    fwhm:
      - 5
      - 10
      - 15
      - 20


We will also add a new variable to point to our input bold file::

    in_bold: '../bids/sub-{subject}/func/sub-{subject}_task-{task}_run-{run}_bold.nii.gz'


In our Snakefile, we then need to replace these hardcoded values with ``config[key]``. The entire updated Snakefile is shown here:

```{literalinclude} step4/Snakefile
  :language: python
```

After these changes, the workflow should still run just like the last step, but now you can make any changes via the config file.

```{asciinema} step4/step4.cast
```


Part II: Snakebids
==================

Now that we have a fully functioning and generic Snakemake workflow, let's see what Snakebids can add.

(step_5)=
Step 5: the bids() function
---------------------------

The first thing we can make use of is the {func}`bids() <snakebids.bids>` function. This provides an easy way to generate bids filenames. This is especially useful when defining output files in your workflow and you have many bids entities.

In our existing workflow, this was our output file:

```py
'results/sub-{subject}/func/sub-{subject}_task-{task}_run-{run}_fwhm-{fwhm}mm_bold.nii.gz'
```

To create the same path using {func}`bids() <snakebids.bids>`, we just need to specify the root directory (`results`), all the bids tags (subject, task, run, fwhm), and the suffix (which includes the extension):

```python
bids(
  root='results',
  subject='{subject}',
  task='{task}',
  run='{run}',
  fwhm='{fwhm}',
  suffix='bold.nii.gz'
)
```


```{note}
To make a snakemake wildcard, we wrapped the `'value'` in curly braces (e.g. `'{value}'`).
```

```{note}
The entities you supply in the {func}`bids() <snakebids.bids>` function do not have to be in the BIDS specification, e.g. ``fwhm`` is not in the spec. But if you do use entities that are in the BIDS specification, snakebids will ensure they go in the correct order (e.g. `sub_*-ses_*-run_*...`). Non-standard entities will be added to the end of the path, just before the suffix, in the order they were defined.
```

The Snakefile with the output filename replaced (in both rules) is below:

```{literalinclude} step5/Snakefile
  :language: python
```


Step 6: parsing the BIDS dataset
--------------------------------

So far, we have had to manually enter the path to input bold file in the config file, and also specify what subjects, tasks, and runs we want processed. Can't we use the fact that we have a BIDS dataset to automate this a bit more?

With Snakemake, there are ways to glob the files to figure out what wildcards are present (e.g. `glob_wildcards`), however, this is not so straightforward with BIDS, since filenames in BIDS often have optional components. E.g. some datasets may have a `ses` tag/sub-directory, and others do not. Also there are often optional user-defined values, such as the `acq` tag, that a workflow in most cases should ignore. Thus, the input that we use in our workflow, `in_bold`, that has wildcards to be generic, would need to be altered for any given BIDS dataset, along with the workflow itself, making this automated BIDS parsing difficult within Snakemake.

Snakebids lets you parse a bids dataset (using [`pybids`](https://bids-standard.github.io/pybids/) under the hood) using a configfile that contain the required wildcards, along with data structures that specify all the wildcard values for all the subjects. This in combination with the {func}`bids() <snakebids.bids>` function, as we will see, can allow one to make snakemake workflows that are compatible with any general bids dataset.

To add this parsing to the workflow, we call the {func}`generate_inputs() <snakebids.generate_inputs()>` function before our rules are defined, and pass along some configuration data to specify the location of the bids directory (`bids_dir`), and the inputs we want to parse for the workflow (`pybids_inputs`). The function returns a {class}`BidsDataset() <snakebids.BidsDataset>`, which we'll assign to a variable called `inputs`:

```py
inputs = snakebids.generate_inputs(
    bids_dir=config['bids_dir'],
    pybids_inputs=config['pybids_inputs'],
    use_bids_lists=True,
)
```

```{note}
Snakebids is transitioning to a new format for {func}`generate_inputs() <snakebids.generate_inputs()>`. Currently, you need to opt-in to the new features by setting `use_bids_lists=True` in `generate_inputs()`, but it will become the default in an upcoming version. We recommend all new users opt-in to maintain long term support, so the tuturial is written using the new syntax. A tutorial for the old syntax can be found on [the v0.5.0 docs](https://snakebids.readthedocs.io/en/v0.5.0/tutorial/tutorial.html#part-ii-snakebids).
```


The config variables we need pre-defined are as follows::

```yaml
bids_dir: '../bids'

pybids_inputs:
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
```

The ``pybids_inputs`` dict defines what types of inputs the workflow can make use of (i.e. the top-level keys, ``bold`` in this case), and for each input, how to filter for them (i.e. the ``filters`` dict), and what BIDS entities to replace with wildcards in the snakemake workflow (i.e. the ``wildcards`` dict).

```{note}
The ``filters`` dict is passed directly to the ``get()`` function in ``pybids``, and thus is quite customizable.
```

%  (TODO: add link)

```{note}
Entries in the `wildcards` list do not have to be in your bids dataset, but if they are, then they will be converted into wildcards (i.e. `task-{task}`) if they are in the filenames. The names for these also correspond with pybids (e.g. `acquisition` maps to `acq`).
```

The {class}`BidsDataset <snakebids.BidsDataset>` class returned by {func}`generate_inputs() <snakebids.generate_inputs>` summarizes the wildcards found in your bids dataset and has parameters to plug those wildcards into your workflow. The first parameter we'll use is {attr}`BidsDataset.path <snakebids.BidsDataset.path>`. You can see what it looks like by printing it after calling {func}`generate_inputs() <snakebids.generate_inputs>`:

```py
print(inputs.path)
```

%TODO: add asciinema
If you ran the workflow, you'd see that `path` is a `dict` with one entry for each named input we specifed under `pybids_config`. In this case, since we only specified one such input: `bold`, `path` has just one item. The values of this dict are generic paths constructed using the wildcards we specified in `pybids_config`. In this case, the path has a `{subject}`, `{session}`, `{acq}`, `{task}`, and `{run}` wildcard. Notice that this path is the same as the path we wrote under `in_bold:` in our `config.yaml` file in {ref}`step 4 <step_4>`. In fact, we can go ahead and replace `config['in_bold']` in our `Snakemake` file with `inputs.path['bold]` and delete `in_bold` from `config.yaml`.

```py
input:
    inputs.path['bold']
```

Step 7: using input wildcards:
------------------------------

{attr}`BidsDataset.path <snakebids.BidsDataset.path>` already grants us a lot of flexibility, but we can still do more! Next we'll look at {attr}`BidsDataset.wildcards <snakebids.BidsDataset.wildcards>`. This is a dict that maps our input names to the brace-wrapped `{wildcards}` we specified in `pybids_config`. If you printed this value in our test workflow, it would look like this:

```py
inputs.wildcards == {
    'bold': {
        'subject': '{subject}',
        'task': '{task}',
        'run': '{run}'
    }
}
```

This is super useful when combined with {func}`bids() <snakebids.bids>`, as we can use the keyword expansion (`**inputs.wildcards[<input_name>]` to set all the wildcard parameters to the {func}`bids() <snakebids.bids>` function. Thus, we can make our workflow even more general, by replacing:

```py
    bids(
        root='results',
        subject='{subject}',
        task='{task}',
        run='{run}',
        fwhm='{fwhm}',
        suffix='bold.nii.gz'
    )
```

with this:

```py
    bids(
        root='results',
        fwhm='{fwhm}',
        suffix='bold.nii.gz',
        **inputs.wildcards['bold']
    )
```


This effectively ensures that any bids entities from the input filenames (that are listed as pybids wildcards) get carried over to the output filenames. Note that we still have the ability to add on additional entities, such as ``fwhm`` here, and set the root directory and suffix.

Another useful config variable is {attr}`entities <snakebids.BidsDataset.entities>`. This is indexed by the pybids input type (e.g. `bold`) and has the list of values that each wildcard can take:

```py
inputs.entities == {
    'bold': {
            'subject': ['001'],
            'task': ['rest'],
            'run': ['1', '2']
    }
}
```

This is a useful dict to use with `expand()` in a target rule, letting us avoid having to specify e.g. the run numbers or task names in the config, and rely on pybids to determine these:

```py
rule all:
    input:
        expand(
            bids(
                root='results',
                fwhm='{fwhm}',
                suffix='bold.nii.gz',
                **config['wildcards']['bold']
            ),
            fwhm=config['fwhm'],
            **inputs.entities['bold']
        )
```


```{note}
Since `expand()` will by default use all combinations of these, this can lead to snakemake searching for inputs that do not exist (e.g. if only one of your subjects has run-2). In this case, you want to use the {attr}`BidsDataset.zip_lists <snakebids.BidsDataset.zip_lists>` instead).
%TODO: maybe just do this in a later step, and/or link to other docs..
```

For reference, here is the updated config file and Snakefile after these changes:

```{literalinclude} step7/config.yml
  :language: yaml
```


```{literalinclude} step7/Snakefile
  :language: python
```

Step 8: creating a command-line executable
------------------------------------------

Now that we have pybids parsing to dynamically configure our workflow inputs based on our BIDS dataset, we are ready to turn our workflow into a [BIDS App](http://bids-apps.neuroimaging.io/), e.g. an app with a standardized command-line interface (e.g. three required positional arguments: ``bids_directory``, ``output_directory``, and ``analysis_level``).

We do this in snakebids by creating an executable python script, which uses the {class}`SnakeBidsApp <snakebids.app.SnakeBidsApp>` class from {mod}`snakebids.app` to run snakemake. An example of this `run.py` script is shown below.

```{literalinclude} step8/run.py
  :language: python
```

However, we will first need to add some additional information to our config file, mainly to define how to parse command-line arguments for this app. This is done with a new `parse_args` dict in the config:

```yaml
parse_args:

  bids_dir:
    help: The directory with the input dataset formatted according
          to the BIDS standard.

  output_dir:
    help: The directory where the output files
          should be stored. If you are running group level analysis
          this folder should be prepopulated with the results of the
          participant level analysis.

  analysis_level:
    help: Level of the analysis that will be performed.
    choices: *analysis_levels

  --participant_label:
    help: The label(s) of the participant(s) that should be analyzed. The label
          corresponds to sub-<participant_label> from the BIDS spec
          (so it does not include "sub-"). If this parameter is not
          provided all subjects should be analyzed. Multiple
          participants can be specified with a space separated list.
    nargs: '+'

  --exclude_participant_label:
    help: The label(s) of the participant(s) that should be excluded. The label
          corresponds to sub-<participant_label> from the BIDS spec
          (so it does not include "sub-"). If this parameter is not
          provided all subjects should be analyzed. Multiple
          participants can be specified with a space separated list.
    nargs: '+'

  --derivatives:
    help: 'Path(s) to a derivatives dataset, for folder(s) that contains multiple derivatives datasets (default: %(default)s) '
    default: False
      type: Path
      nargs: '+'
```


The above is standard boilerplate for any BIDS app.  You can also define any new command-line arguments you wish. Snakebids uses the {mod}`argparse` module, and each entry in this `parse_args` dict thus becomes a call to {meth}`add_argument() <argparse.ArgumentParser.add_argument>` from {class}`argparse.ArgumentParser`. When you run the workflow, snakebids adds the named argument values to the config dict, so your workflow can make use of it as if you had manually added the variable to your configfile.

Arguments that will receive paths should be given the item `type: Path`, as is done for `--derivatives` in the example above. Without this annotation, paths given to keyword arguments will be interpreted relative to the output directory. Indicating `type: Path` will tell Snakebids to first resolve the path according to your current working directory.

```{note}
`bids_dir` and `output_dir` are always resolved relative to the current directory and do not need any extra annotation.
```

BIDS apps also have a requried `analysis_level` positional argument, so there are some config variables to set this as well. The analysis levels are in an `analysis_levels` list in the config, and also as keys in a `targets_by_analysis_level` dict, which can be used to map each analysis level to the name of a target rule::

```yaml
targets_by_analysis_level:
  participant:
    - ''  # if '', then the first rule is run

analysis_levels: &analysis_levels
  - participant
```

Note: since we specified a `''` for the target rule, no target rule will be specified, so snakemake will just default to the first rule in the workflow.

Make the `run.py` script executable (`chmod a+x run.py`) and try running it now.

%TODO: add asciinema

The updated configfile is here (the Snakefile did not change in this step):

```{literalinclude} step8/config.yml
  :language: yaml
```
