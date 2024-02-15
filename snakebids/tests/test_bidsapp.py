# ruff: noqa: PLR2004
from __future__ import annotations

import argparse
from typing import Any

import pytest

from snakebids import bidsapp, plugins
from snakebids.bidsapp.args import ArgumentGroups


class TestRunner:
    def setup_method(self):
        self.hooks_run = 0

    @bidsapp.hookimpl
    def initialize_config(self, config: dict[str, Any]):
        self.hooks_run += 1
        assert not config
        config["initialized"] = True

    @bidsapp.hookimpl
    def add_cli_arguments(
        self,
        parser: argparse.ArgumentParser,
        config: dict[str, Any],
        argument_groups: ArgumentGroups,
    ):
        self.hooks_run += 1
        assert len(config) == 1
        assert config["initialized"]
        assert len(parser._actions) == 1
        assert len(argument_groups) == 2
        config["added_args"] = True
        parser.add_argument("arg_one")
        parser.add_argument("--arg-two")

    @bidsapp.hookimpl
    def handle_unknown_args(self, args: list[str], config: dict[str, Any]):
        self.hooks_run += 1
        assert args == ["unknown"]
        assert len(config) == 2
        assert config["initialized"]
        assert config["added_args"]
        config["unknown_args"] = args

    @bidsapp.hookimpl
    def update_cli_namespace(self, namespace: dict[str, Any], config: dict[str, Any]):
        self.hooks_run += 1
        assert len(config) == 3
        assert config["initialized"]
        assert config["added_args"]
        assert config["unknown_args"] == ["unknown"]
        assert len(namespace) == 2
        assert namespace["arg_one"] == "known"
        assert namespace["arg_two"] == "known"
        config["arg_two"] = namespace.pop("arg_two")

    @bidsapp.hookimpl
    def finalize_config(self, config: dict[str, Any]):
        self.hooks_run += 1
        assert len(config) == 5
        assert config["initialized"]
        assert config["added_args"]
        assert config["unknown_args"] == ["unknown"]
        assert config["arg_one"] == "known"
        assert config["arg_two"] == "known"
        config["finalized"] = True

    @bidsapp.hookimpl
    def run(self, config: dict[str, Any]):
        self.hooks_run += 1
        assert len(config) == 6
        assert config["initialized"]
        assert config["added_args"]
        assert config["unknown_args"] == ["unknown"]
        assert config["arg_one"] == "known"
        assert config["arg_two"] == "known"
        assert config["finalized"]
        config["run"] = True

    def test_build_parser(self):
        app = bidsapp.app(plugins=[self])
        app.build_parser()
        assert self.hooks_run == 2
        app.build_parser()
        assert self.hooks_run == 2

    def test_parse_args(self):
        app = bidsapp.app(plugins=[self])
        app.parse_args(args=["known", "--arg-two", "known", "unknown"])
        assert self.hooks_run == 5
        app.parse_args(args=["known", "--arg-two", "known", "unknown"])
        assert self.hooks_run == 5

    def test_run(self):
        app = bidsapp.app(plugins=[self])
        app.run(args=["known", "--arg-two", "known", "unknown"])
        assert self.hooks_run == 6
        app.run(args=["known", "--arg-two", "known", "unknown"])
        assert self.hooks_run == 6


class TestDependencies:
    DEPENDENCIES = (plugins.BidsArgs(),)

    def setup_method(self):
        self.called = False

    @bidsapp.hookimpl
    def add_cli_arguments(
        self,
        parser: argparse.ArgumentParser,
    ):
        self.called = True
        assert len(parser._actions) > 1

    def test_dependencies_are_added(self):
        app = bidsapp.app(plugins=[self])
        assert len(app.pm.list_name_plugin()) == 2

    def test_dependences_registered_after(self):
        app = bidsapp.app(plugins=[self])
        app.build_parser()
        assert self.called

    def test_dependencies_wont_register_twice(self):
        app = bidsapp.app(
            plugins=[
                self,
                plugins.BidsArgs(
                    bids_dir=False,
                    output_dir=False,
                    analysis_level=False,
                    participant_label=False,
                    exclude_participant_label=False,
                    # leave derivatives in so that the DependenciesPlugin
                    # assertion doesn't get tripped
                ),
            ]
        )
        app.run()
        assert self.called
        assert len(app.parser._actions) == 2


class TestDependencyDeduplication:
    DEPENDENCIES = (
        plugins.BidsArgs(),
        plugins.ComponentEdit(),
        plugins.CliConfig(),
        plugins.Pybidsdb(),
        plugins.BidsValidator(),
    )

    def get_plugins(self):
        # We need to regenerate all instances again to check that the `isinstance`
        # equality checks in all plugins are working
        return (
            plugins.BidsArgs(),
            plugins.ComponentEdit(),
            plugins.CliConfig(),
            plugins.Pybidsdb(),
            plugins.BidsValidator(),
        )

    def test_duplicate_normally_caught(self):
        with pytest.raises(ValueError, match="Plugin already registered"):
            bidsapp.app(plugins=[*self.get_plugins(), *self.DEPENDENCIES])

    def test_no_duplicate_plugins_added_when_dependencies(self):
        plugins = self.get_plugins()
        app = bidsapp.app(plugins=[self, *plugins])
        assert len(app.pm.list_name_plugin()) == len(plugins) + 1


class TestArgumentGroups:
    def setup_method(self):
        self.group_names = None

    @bidsapp.hookimpl
    def add_cli_arguments(
        self,
        argument_groups: ArgumentGroups,
    ):
        assert argument_groups["new"]
        self.group_names = list(argument_groups)

    def test_argument_groups_visible_to_add_arg_hook(self):
        app = bidsapp.app(plugins=[self])
        app.parser.add_argument_group("foo")
        app.run()
        assert self.group_names
        assert set(self.group_names).issuperset({"foo", "new"})
        assert app.argument_groups["new"].title == "new"
