from __future__ import annotations

import textwrap
import warnings

from colorama import Fore, Style

WARN_TEMPLATE = f"""\
{Fore.YELLOW}{Style.BRIGHT}[{{category}}]{Style.NORMAL} {{filename}}:{{lineno}} \
{Fore.RESET}{{line}}
{{message}}

"""


def formatwarning(
    message: Warning | str,
    category: type[Warning],
    filename: str,
    lineno: int,
    line: str | None,
):
    """Format warning messages."""
    if line is None:
        try:
            import linecache

            line = linecache.getline(filename, lineno).rstrip()
        except Exception:  # pragma: no cover # noqa: BLE001
            # When a warning is logged during Python shutdown, linecache
            # and the import machinery don't work anymore
            line = None
            linecache = None

    return WARN_TEMPLATE.format(
        message=textwrap.indent(
            message.args[0] if isinstance(message, Warning) else message,
            prefix="  ",
        ),
        category=category.__name__,
        filename=filename,
        lineno=lineno,
        line=line,
    )


warnings.formatwarning = formatwarning
