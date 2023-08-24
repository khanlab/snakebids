from __future__ import annotations

import re
import subprocess
import sys
from pathlib import Path

import jinja2.parser
from jinja2 import nodes
from jinja2.ext import Extension


class GitConfigExtension(Extension):
    tags = {"gitconfig"}  # noqa: RUF012
    _config: dict[str, str]

    def __init__(self, env: jinja2.Environment) -> None:
        self._config = {}

        try:
            config_list = subprocess.check_output(
                [executable(), "config", "-l"], stderr=subprocess.STDOUT
            ).decode()

            m = re.findall("(?ms)^([^=]+)=(.*?)$", config_list)
            if m:
                for group in m:
                    self._config[group[0]] = group[1]
        except (subprocess.CalledProcessError, OSError):
            pass

    def get(self, key: str, default: str | None = None) -> str | None:
        return self._config.get(key, default)

    def __getitem__(self, item: str) -> str:
        return self._config[item]

    def parse(self, parser: jinja2.parser.Parser):
        lineno = next(parser.stream).lineno

        node = parser.parse_expression()

        if not isinstance(node, nodes.Const):
            raise ValueError("Argument to `gitconfig` must be a string")
        call_method = self.call_method(
            "get",
            [node],
            lineno=lineno,
        )
        return nodes.Output([call_method], lineno=lineno)


def executable() -> str:
    _executable = None

    if sys.platform == "win32":
        # Finding git via where.exe
        where = "%WINDIR%\\System32\\where.exe"
        paths = subprocess.check_output(
            [where, "git"], shell=True, encoding="oem"
        ).split("\n")
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
        raise RuntimeError("Unable to find a valid git executable")

    return _executable
