import json

import jinja2
from jinja2.ext import Extension


def toml_string(item: str):
    return json.dumps(item, ensure_ascii=False).replace("\x7F", "\\u007f")


class TomlEncodeExtension(Extension):
    def __init__(self, env: jinja2.Environment):
        env.filters["toml_string"] = toml_string  # type: ignore
