# Plugins

Plugins are a feature in Snakebids that enables you to add custom functionality to your Snakebids application after parsing CLI arguments but before invoking Snakemake. For example, you can use a plugin to perform BIDS validation of your Snakebids app's input, which ensures your app is only executed if the input dataset is valid. You can either use those that are distributed with Snakebids (see [Using plugins](#using-plugins)) or create your own plugins (see [Creating plugins](#creating-plugins)).

## Using plugins
To add one or more plugins to your {class}`SnakeBidsApp <snakebids.app.SnakeBidsApp>`, pass them to the {class}`~snakebids.app.SnakeBidsApp` constructor via the {attr}`~snakebids.app.SnakeBidsApp.plugins` parameter. Your plugin will have access to CLI parameters (after they've been parsed) via their names in {attr}`SnakeBidsApp.config <snakebids.app.SnakeBidsApp.config>`. Any modifications to that config dictionary made by the plugin will be carried forward into the workflow.

As an example, the distributed {class}`BidsValidator` plugin can be used to run the [BIDS Validator](https://github.com/bids-standard/bids-validator) on the input directory like so:

```py
import subprocess

from snakebids.app import SnakeBidsApp
from snakebids.plugins.validator import BidsValidator

SnakeBidsApp(
    "path/to/snakebids/app",
    plugins=[BidsValidator]
).run_snakemake()
```

## Creating plugins
A plugin is a function or callable class that accepts a {class}`SnakeBidsApp <snakebids.app.SnakeBidsApp>` as input and returns a modified {class}`SnakeBidsApp` or `None`.

As an example, a simplified version of the bids-validator plugin that runs the [BIDS Validator](https://github.com/bids-standard/bids-validator) could be defined as follows:

```py
import subprocess 

from snakebids.app import SnakeBidsApp
from snakebids.exceptions import SnakeBidsPluginError

class BidsValidator:
    """Perform BIDS validation of dataset

    Parameters
    -----------
    app
        Snakebids application to be run
    """

    def __call__(self, app: SnakeBidsApp) -> None:
        # Skip bids validation
        if app.config["plugins.validator.skip"]:
            return

        try:
            subprocess.run(["bids-validator", app.config["bids_dir"]], check=True)
        except subprocess.CalledProcessError as err:
            raise InvalidBidsError from err


class InvalidBidsError(SnakebidsPluginError):
    """Error raised if input BIDS dataset is invalid, 
    inheriting from SnakebidsPluginError.
    """
```

```{note}
When adding plugin-specific parameters to the config dictionary, it is recommended to use namespaced keys (e.g. ``plugins.validator.skip``). This will help ensure plugin-specific parameters do not conflict with other parameters already defined in the dictionary or by other plugins.
```


A plugin can be used to implement any logic that can be handled by a Python function. In the above example, you may also want to add some logic to check if the BIDS Validator is installed and pass along a custom error message if it is not. Created plugins can then be used within a Snakebids workflow, similar to the example provided in [Using plugins](#using-plugins) section. Prospective plugin developers can take a look at the source of the `~snakebids.plugins` module for examples.


```{note}
When creating a custom error for your Snakebids plugin, it is recommended to inherit from {class}`SnakebidsPluginError <snakebids.exceptions.SnakebidsPluginError>` such that errors will be recognized as a plugin error.
```