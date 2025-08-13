from __future__ import annotations

import jinja2.parser
from colorama import Fore
from jinja2.ext import Extension


class ColoramaExtension(Extension):
    """Include colorama foreground colors in the global environment."""

    def __init__(self, env: jinja2.Environment):
        super().__init__(env)
        env.globals["Fore"] = Fore  # type: ignore
