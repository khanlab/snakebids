# Overview

Snakebids apps rely on a configuration file (`snakebids.yml`). This file specifies which files from a BIDS dataset should be used as input. The apps also utilize workflow definitions, which are written in one or more Snakefile(s) and determine how the input files are processed.

```{note}
For an easy setup of new Snakebids apps with convenient command-line functions, we recommend installing Snakebids using `pipx`. Visit the [following page](https://pypa.github.io/pipx/) for instructions on how to install `pipx`.
```

Once Snakebids is installed, you can generate a customized Snakebids project by running the command `snakebids create` and providing the necessary information when prompted.
