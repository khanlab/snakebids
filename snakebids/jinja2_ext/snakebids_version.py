from __future__ import annotations

import importlib.metadata as impm
import re

import jinja2.parser
import requests
from jinja2.ext import Extension


class SnakebidsVersionExtension(Extension):
    """Retrieve the latest snakebids vesion from pypi.

    Stores value in the ``snakebids_version`` global variable.
    """

    def __init__(self, env: jinja2.Environment):
        env.globals["snakebids_version"] = self._lookup_version()  # type: ignore
        super().__init__(env)

    def _lookup_version(self):
        request = requests.get("https://pypi.org/pypi/snakebids/json", timeout=10)
        version_regex = re.compile(r"\d+\.\d+\.\d+")
        try:
            request.raise_for_status()
            version = request.json()["info"]["version"]
            version = version_regex.fullmatch(version).group()  # type: ignore
        except (requests.HTTPError, KeyError, AttributeError, TypeError):
            version = impm.version("snakebids")
            if version == "0.0.0" or not version_regex.fullmatch(version):
                return None
        return version
