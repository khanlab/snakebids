from __future__ import annotations

import importlib.metadata as impm
import re

import jinja2.parser
import requests
from jinja2.ext import Extension


class SnakebidsVersionExtension(Extension):
    def __init__(self, env: jinja2.Environment):
        env.globals["snakebids_version"] = self._lookup_version()  # type: ignore
        super().__init__(env)

    def _lookup_version(self):
        request = requests.get("https://pypi.org/pypi/snakebids/json")
        version_regex = re.compile(r"\d+\.\d+\.\d+")
        try:
            request.raise_for_status()
            version = request.json()["info"]["version"]
            if not re.fullmatch(version_regex, version):
                raise TypeError()
            return version
        except (requests.HTTPError, KeyError, TypeError):
            version = impm.version("snakebids")
            if version == "0.0.0" or not re.fullmatch(version_regex, version):
                return None
            return version
