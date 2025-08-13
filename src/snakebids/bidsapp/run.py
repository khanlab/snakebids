from __future__ import annotations

import argparse
import sys
from collections import defaultdict
from typing import TYPE_CHECKING, Any, Iterable, Protocol, Sequence, TypedDict

import attrs
import pluggy
from attrs.converters import default_if_none
from typing_extensions import TypeAlias, Unpack

from snakebids.bidsapp import hookspecs
from snakebids.bidsapp.args import ArgumentGroups

if TYPE_CHECKING:
    from argparse import _ArgumentGroup, _FormatterClass

    from pluggy._hooks import _Plugin

if sys.version_info >= (3, 9):

    class ArgumentParserArgs(TypedDict, total=False):
        """Arguments passed on to :class:`argparse.ArgumentParser`."""

        prog: str | None
        usage: str | None
        description: str | None
        epilog: str | None
        parents: Sequence[argparse.ArgumentParser]
        formatter_class: _FormatterClass
        prefix_chars: str
        fromfile_prefix_chars: str | None
        argument_default: Any
        conflict_handler: str
        add_help: bool
        allow_abbrev: bool
        exit_on_error: bool
else:

    class ArgumentParserArgs(TypedDict, total=False):
        """Arguments passed on to :class:`argparse.ArgumentParser`."""

        prog: str | None
        usage: str | None
        description: str | None
        epilog: str | None
        parents: Sequence[argparse.ArgumentParser]
        formatter_class: _FormatterClass
        prefix_chars: str
        fromfile_prefix_chars: str | None
        argument_default: Any
        conflict_handler: str
        add_help: bool
        allow_abbrev: bool


def app(
    plugins: Iterable[_Plugin] | None = None,
    *,
    config: dict[str, Any] | None = None,
    **argparse_args: Unpack[ArgumentParserArgs],
) -> _Runner:
    """Create a BIDSApp.

    Parameters
    ----------
    plugins
        List of snakebids plugins to apply to app (see :ref:`using-plugins`).
    config
        Initial config. Will be updated with parsed arguments from ``parser`` and
        potentially modified by plugins
    **argparse_args
        Arguments passed on transparently to :class:`argparse.ArgumentParser`.
    """
    pm = pluggy.PluginManager("snakebids")
    pm.add_hookspecs(hookspecs)
    parser = argparse.ArgumentParser(**argparse_args)
    if plugins is not None:
        _add_plugin_dependencies(pm, plugins, registered_ok=False)
    return _Runner(pm=pm, parser=parser, config=config)


class _DependingPlugin(Protocol):
    DEPENDENCIES: tuple[_Plugin, ...]


def _add_plugin_dependencies(
    pm: pluggy.PluginManager,
    plugins: Iterable[_Plugin | _DependingPlugin],
    *,
    registered_ok: bool,
) -> None:
    dependencies: list[_Plugin] = []
    for plugin in plugins:
        if registered_ok and pm.is_registered(plugin):
            continue
        pm.register(plugin)
        if hasattr(plugin, "DEPENDENCIES"):
            dependencies.extend(plugin.DEPENDENCIES)  # type: ignore
    if not dependencies:
        return None
    return _add_plugin_dependencies(pm, dependencies, registered_ok=True)


if TYPE_CHECKING:
    _ArgumentGroupDictType: TypeAlias = "defaultdict[str, _ArgumentGroup]"
else:
    _ArgumentGroupDictType = defaultdict


class _ArgumentGroupDict(_ArgumentGroupDictType):
    """Dict that converts missing keys into argument groups."""

    def __init__(self, parser: argparse.ArgumentParser):
        self.parser = parser

    def __missing__(self, key: str):
        new = self.parser.add_argument_group(key)
        self[key] = new
        return new


@attrs.define
class _Runner:
    """Runtime manager for BIDS apps.

    Manages the parser, config, and plugin relay for snakebids BIDS apps. This class
    should not be constructed directly, but built using the
    :func:`bidsapp.app() <snakebids.bidsapp.app>` function.
    """

    pm: pluggy.PluginManager
    """Reference to the plugin manager."""

    parser: argparse.ArgumentParser
    """Parser used for parsing CLI arguments.

    The parser may be manipulated before running action methods to add initial arguments
    (via :meth:`parser.add_argument() <argparse.ArgumentParser.add_argument>`) and
    argument groups (via :meth:`parser.add_argument_group()
    <argparse.ArgumentParser.add_argument_group>`). Any argument groups added will
    automatically be indexed in :attr:`.argument_groups` and made available to plugins
    via their ``title``.
    """

    config: dict[str, Any] = attrs.field(converter=default_if_none(factory=dict))
    """Configuration dictionary for passing data between plugins."""

    argument_groups: ArgumentGroups = attrs.field(
        default=attrs.Factory(
            lambda self: _ArgumentGroupDict(self.parser), takes_self=True
        )
    )
    """Argument group reference accessible to plugins for organizing CLI arguments.

    Any groups added to the :attr:`parser` before the action methods are called will
    be automatically added, indexed by their titles. Thus, this object should typically
    only be used within plugins.
    """

    _parser_built: bool = attrs.field(default=False, init=False)
    _processed: bool = attrs.field(default=False, init=False)
    _run: bool = attrs.field(default=False, init=False)

    def build_parser(self):
        """Run plugins affecting the ``parser`` without yet parsing arguments."""
        if not self._parser_built:
            self.pm.hook.initialize_config(config=self.config)
            for group in self.parser._action_groups:  # noqa: SLF001
                if group.title is not None:
                    self.argument_groups[group.title] = group
            self.pm.hook.add_cli_arguments(
                parser=self.parser,
                argument_groups=self.argument_groups,
                config=self.config,
            )
            self._parser_built = True
        return self

    def parse_args(self, args: list[str] | None = None):
        """Run all plugins and parse arguments."""
        if not self._processed:
            self.build_parser()
            args = sys.argv[1:] if args is None else args
            argv: list[str] | None = self.pm.hook.get_argv(
                argv=args, config=self.config
            )
            namespace, unknown = self.parser.parse_known_args(
                args=args if argv is None else argv
            )
            self.pm.hook.handle_unknown_args(args=unknown, config=self.config)
            self.pm.hook.update_cli_namespace(
                namespace=namespace.__dict__, config=self.config
            )
            self.config.update(namespace.__dict__)
            self.pm.hook.finalize_config(config=self.config)
            self._processed = True
        return self

    def run(self, args: list[str] | None = None):
        if not self._run:
            self.parse_args(args)
            self.pm.hook.run(config=self.config)
            self._run = True
