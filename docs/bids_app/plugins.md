# Plugins

Plugins are a Snakebids feature that allow you to add arbitrary behaviour to your Snakebids app after CLI arguments are parsed but before Snakemake is invoked. For example, you might add BIDS validation of an input dataset to your app via a plugin, so your app is only run if the input dataset is valid.

A plugin is simply a function that takes a {class}`SnakeBidsApp <snakebids.app.SnakeBidsApp>` as input and returns either a modified {class}`SnakeBidsApp <snakebids.app.SnakeBidsApp>` or `None`. To add one or more plugins to your {class}`SnakeBidsApp <snakebids.app.SnakeBidsApp>`, use the method {func}`add_plugins <snakebids.app.SnakeBidsApp.add_plugins>`, which will typically be easiest in a `run.py` or similar entrypoint to your app. Your plugin will have access to CLI parameters via `SnakeBidsApp.config`, but you're responsible for making sure that your app provides the properties a plugin expects to find in the config.
