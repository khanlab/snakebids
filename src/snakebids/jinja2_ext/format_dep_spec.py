import jinja2
from jinja2.ext import Extension

from snakebids.jinja2_ext.toml_encode import toml_string


def format_poetry(item: str):
    """Format a pip style dependency specification for a poetry pyproject.toml.

    Only supports urls (prefixed with @) and version specifiers. No markers or extras.
    The package name should already be stripped.
    """
    if item.strip().startswith("@"):
        url = item[1:].strip()
        if url.startswith("git+"):
            return f"{{ git = {toml_string(url[4:])} }}"
        if url.startswith("file://"):
            return f"{{ path = {toml_string(url[7:])} }}"
        return f"{{ url = {toml_string(url)} }}"
    return toml_string(item)


class FormatDepSpec(Extension):
    """Enable the toml_string filter, which encodes provided value as toml."""

    def __init__(self, env: jinja2.Environment):
        env.filters["format_poetry"] = format_poetry  # type: ignore
