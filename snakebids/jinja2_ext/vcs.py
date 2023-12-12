from __future__ import annotations

import re
import subprocess as sp
import sys
from pathlib import Path

import jinja2.parser
from jinja2 import nodes
from jinja2.ext import Extension
from typing_extensions import override


class GitConfigExtension(Extension):
    """Retrieve settings from git configuration.

    Shortcode activated with ``gitconfig`` plus the key name::

        {% gitconfig "user.name" %}
    """

    tags = {"gitconfig"}  # noqa: RUF012
    _config: dict[str, str]

    def __init__(self, env: jinja2.Environment) -> None:
        self._config = {}

        try:
            config_list = sp.check_output(
                [executable(), "config", "-l"], stderr=sp.STDOUT
            ).decode()

            m = re.findall("(?ms)^([^=]+)=(.*?)$", config_list)
            if m:
                for group in m:
                    self._config[group[0]] = group[1]
        except (sp.CalledProcessError, OSError):
            pass

    def get(self, key: str, default: str | None = None) -> str | None:
        """Get specified configuration key."""
        return self._config.get(key, default)

    def __getitem__(self, item: str) -> str:
        return self._config[item]

    @override
    def parse(self, parser: jinja2.parser.Parser):
        lineno = next(parser.stream).lineno

        node = parser.parse_expression()

        if not isinstance(node, nodes.Const):
            msg = "Argument to `gitconfig` must be a string"
            raise TypeError(msg)
        call_method = self.call_method(
            "get",
            [node],
            lineno=lineno,
        )
        return nodes.Output([call_method], lineno=lineno)


def executable() -> str:
    """Retrieve os-dependent command for git."""
    _executable = None

    if sys.platform == "win32":
        # Finding git via where.exe
        where = "%WINDIR%\\System32\\where.exe"
        paths = sp.check_output([where, "git"], encoding="oem").split("\n")
        for path in paths:
            if not path:
                continue

            _path = Path(path.strip())
            try:
                _path.relative_to(Path.cwd())
            except ValueError:
                _executable = str(_path)

                break
    else:
        _executable = "git"

    if _executable is None:  # type: ignore
        msg = "Unable to find a valid git executable"
        raise RuntimeError(msg)

    return _executable
