import json

import jinja2
from jinja2.ext import Extension


def toml_string(item: str):
    """Encode string for inclusion in toml.

    Technically encodes as json, a (mostly) strict subset of toml, with some encoding
    fixes
    """
    return json.dumps(item, ensure_ascii=False).replace("\x7F", "\\u007f")


class TomlEncodeExtension(Extension):
    """Enable the toml_string filter, which encodes provided value as toml."""

    def __init__(self, env: jinja2.Environment):
        env.filters["toml_string"] = toml_string  # type: ignore
