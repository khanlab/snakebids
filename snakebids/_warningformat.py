from __future__ import annotations

import textwrap
import warnings
from pathlib import Path

from colorama import Fore, Style

WARN_TEMPLATE = f"""\
{Fore.RED}{Style.BRIGHT}[{{category}}]{Style.NORMAL} {{filename}}:{{lineno}} \
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
        with Path(filename).open() as f:
            for _ in range(lineno):
                line = f.readline().strip()

    return WARN_TEMPLATE.format(
        message=textwrap.indent(
            "\n".join(
                "\n".join(
                    textwrap.wrap(
                        line,
                        width=80,
                    )
                )
                for line in (
                    message.args[0] if isinstance(message, Warning) else message
                ).splitlines()
            ),
            prefix="  ",
        ),
        category=category.__name__,
        filename=filename,
        lineno=lineno,
        line=line,
    )


warnings.formatwarning = formatwarning
