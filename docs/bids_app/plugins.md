# Plugins

Plugins allow you to extend the functionality of your [BIDS app](#snakebids.bidsapp) before and after parsing CLI arguments. For example, you can use a plugin to perform BIDS validation of your Snakebids app's input, which ensures your app is only executed if the input dataset is valid. You can either use those that are distributed with Snakebids (see [Using plugins](#using-plugins)) or create your own plugins (see [Creating plugins](#creating-plugins)).

```{note}
For a full list of plugins distributed with Snakebids, see the [Plugins reference](/api/plugins) page.
```

Nearly all of the functionality provided by {mod}`snakebids.bidsapp` is provided by plugins, including {class}`~snakebids.plugins.snakemake.SnakemakeBidsApp`.

Unlike in libraries such as pytest, plugins must be explicitly enabled to be included in the app. Installing them with pip is not enough! This affords a great deal of control over how the BIDS app is executed.

(using-plugins)=
## Using plugins
To add plugins to your {class}`bidsapp <snakebids.app.SnakeBidsApp>`, pass them to the {class}`~snakebids.app.SnakeBidsApp` constructor via the {attr}`~snakebids.app.SnakeBidsApp.plugins` parameter. Plugins are executed in LIFO order (last in, first out).

As an example, the {class}`~snakebids.plugins.BidsValidator` plugin can be used to run the [BIDS Validator](https://github.com/bids-standard/bids-validator) on the input directory like so:

```python
from snakebids import bidsapp, plugins

bidsapp.app([
    plugins.SnakemakeBidsApp("path/to/snakebids/app"),
    plugins.BidsValidator,
]).run()
```

(dependencies)=
### Dependencies

Some plugins depend on other plugins. These dependencies will be loaded even if not specified in {func}`bidsapp.app <snakebids.bidsapp.app>`. If a dependency needs explicit configuration, however, they may still be safely provided in app initialization, as snakebids will prevent duplicate registration.

```{currentmodule} snakebids.plugins
```

For example, {class}`SnakemakeBidsApp` depends on {class}`BidsArgs`, {class}`CliConfig`, {class}`ComponentEdit`, etc, but these plugins may still be specified to change their configuration. The order of specification will be respected (e.g. [LIFO](inv:#callorder)).

```python
from snakebids import bidsapp, plugins

bidsapp.app([
    plugins.SnakemakeBidsApp("path/to/snakebids/app"),
    # specify the BidsArgs plugin to override the argument group
    plugins.BidsArgs(argument_group="MAIN"),
]).run()
```

(creating-plugins)=
## Creating plugins

Plugins are implemented using [`pluggy`](inv:pluggy#index), the plugin system used and maintained by pytest. Actions are executed in one of several hooks. These are called at specified times as the argument parser is built, arguments are parsed, and the config is formatted.

A plugin is a class or module with methods or functions wrapped with the snakebids plugin hook decorator: {attr}`snakebids.bidsapp.hookimpl`. The name of the function determines the stage of app initialization at which it will be called. Each recognized function name (known as specs) comes with a specified set of available arguments. Not all the available arguments need to be used, however, they must be given the correct name. The [API documentation](#snakebids.bidsapp.hookimpl) contains the complete list of available specs, their corresponding initialization stages, and the arguments they can access.

As an example, a simplified version of the bids-validator plugin that runs the [BIDS Validator](https://github.com/bids-standard/bids-validator) could be defined as follows:

```python
import argparse
import subprocess
from typing import Any

from snakebids import bidsapp

class BidsValidator:
    """Perform BIDS validation of dataset

    Parameters
    -----------
    app
        Snakebids application to be run
    """

    @bidsapp.hookimpl
    def add_cli_arguments(self, parser: argparse.ArgumentParser):
        parser.add_argument("--skip-validation", dest="plugins.validator.skip")

    @bidsapp.hookimpl
    def finalize_config(self, config: dict[str, Any]) -> None:
        # Skip bids validation
        if config["plugins.validator.skip"]:
            return

        try:
            subprocess.run(
                ["bids-validator", config["bids_dir"]], check=True
            )
        except subprocess.CalledProcessError as err:
            raise InvalidBidsError from err


class InvalidBidsError(SnakebidsPluginError):
    """Error raised if input BIDS dataset is invalid,
    inheriting from SnakebidsPluginError.
    """
```

```{currentmodule} snakebids.bidsapp.hookspecs
```

In this example, two hooks were used. The {func}`add_cli_arguments` hook is called before CLI arguments are parsed. Here, it adds an argument allowing end users of our app to skip bids validation. Note that the [spec](#snakebids.bidsapp.hookspecs.add_cli_arguments) specifies three arguments available to this hook (`parser`, `config`, and `argument_groups`), however, we only used `parser` here.

```{note}
When adding plugin-specific parameters to the config dictionary, it is recommended to use namespaced keys (e.g. ``plugins.validator.skip``). This will help ensure plugin-specific parameters do not conflict with other parameters already defined in the dictionary or by other plugins.
```

The {func}`finalize_config` hook is called after `config` is updated with the results of argument parsing. In our plugin, this is where validation is actually performed. Note how the argument added in the previous hook is now read from `config`. Any modifications made to `config` will be carried forward into the app's remaining lifetime.

A plugin can be used to implement any logic that can be handled by a Python function. In the above example, you may also want to add some logic to check if the BIDS Validator is installed and pass along a custom error message if it is not. Created plugins can then be used within a Snakebids workflow, similar to the example provided in [Using plugins](#using-plugins) section. Prospective plugin developers can take a look at the source of the `snakebids.plugins` module for examples.

```{note}
When creating a custom error for your Snakebids plugin, it is recommended to inherit from {class}`SnakebidsPluginError <snakebids.exceptions.SnakebidsPluginError>` such that errors will be recognized as a plugin error.
```

### Specifying dependencies

[Dependencies](#dependencies) may be specified using a `DEPENDENCIES` attribute in your plugin class or module. It should be set to a tuple containing fully initialized plugin references. These dependencies will be registered after the depending plugin and therefore run first (due to pluggy's [LIFO order](inv:#callorder)).

```python
from snakebids import bidsapp, plugins

class MyPlugin:
    DEPENDENCIES = (
        plugins.CliConfig(),
        plugins.BidsArgs(),
    )

    @bidsapp.hookimpl
    def finalize_config(self, config):
        ...
```
