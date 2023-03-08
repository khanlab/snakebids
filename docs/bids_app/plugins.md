# Plugins

Plugins are a Snakebids feature that allow you to add arbitrary behaviour to your Snakebids app after CLI arguments are parsed but before Snakemake is invoked. For example, you might add BIDS validation of an input dataset to your app via a plugin, so your app is only run if the input dataset is valid.

A plugin is simply a function that takes a {class}`SnakeBidsApp <snakebids.app.SnakeBidsApp>` as input and returns either a modified {class}`SnakeBidsApp <snakebids.app.SnakeBidsApp>` or `None`. To add one or more plugins to your {class}`SnakeBidsApp <snakebids.app.SnakeBidsApp>`, use the method {func}`add_plugins <snakebids.app.SnakeBidsApp.add_plugins>`, which will typically be easiest in a `run.py` or similar entrypoint to your app. Your plugin will have access to CLI parameters (after they've been parsed) via their names in {attr}`SnakeBidsApp.config <snakebids.app.SnakeBidsApp.config>`. Any modifications to that config dictionary will be carried forward into the workflow.

As an example, a plugin could run the [BIDS Validator](https://github.com/bids-standard/bids-validator) on the input directory like so:

```
import subprocess

from snakebids.app import SnakeBidsApp

def bids_validate(app: SnakeBidsApp) -> None:
    if app.config["skip_bids_validation"]:
        return

    try:
        subprocess.run(["bids-validator", app.config["bids_dir"]], check=True)
    except subprocess.CalledProcessError as err:
        raise InvalidBidsError from err

class InvalidBidsError(Exception):
    """Error raised if an input BIDS dataset is invalid."""

app = SnakeBidsApp("path/to/snakebids/app")
app.add_plugins([bids_validate])
app.run_snakemake()

```

You would also want to add some logic to check if the BIDS Validator is installed and pass along the error message, but the point is that a plugin can do anything that can be handled by a Python function.
