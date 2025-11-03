"""Recompile the bids function stub file based on latest specs."""

from __future__ import annotations

import inspect
import itertools as it
import re
from collections.abc import Iterable
from pathlib import Path
from types import ModuleType

from snakebids.paths import _config, specs
from snakebids.paths._templates import spec_func
from snakebids.paths._utils import BidsPathSpecFile, get_specs


def generate_stub(mod: ModuleType, imports: list[str], funcs: Iterable[str]):
    """Write stub file to provided module.

    Parameters
    ----------
    mod
        module for which stub file should be written
    imports
        list of imports and other statements to appear at the beginning of the file
    funcs
        list of function declarations
    """
    prelude = [
        "# This stub file is automatically generated",
        "# It can be updated using::",
        "#",
        "#      poetry run poe update_bids",
        "",
    ]
    if mod.__file__ is None:
        msg = f"{mod} has no associated file"
        raise FileNotFoundError(msg)
    Path(mod.__file__).with_suffix(".pyi").write_text(
        "\n".join(it.chain(prelude, imports, funcs))
    )


SOURCE_TEMPLATE = """
# <AUTOUPDATE>
# The code between these tags is automatically generated. Do not
# manually edit
# To update, run::
#
#       poetry run poe update-bids
#
{all_block}
{statements}
# </AUTOUPDATE>
"""

ALL_TEMPLATE = """
if not TYPE_CHECKING:
    __all__ = [ # noqa:F822
        {all}
    ]

    def __dir__():
        return __all__
"""


def update_source(
    mod: ModuleType,
    all_items: Iterable[str] | None,
    statements: list[str] | None = None,
):
    """Update python source file with __all__ declaration and provided statements.

    Parameters
    ----------
    mod
        Module to update
    all_items
        Items to include in __all__ declaration
    statements
        Statements to include after __all__ declaration block
    """
    source = inspect.getsource(mod)
    compiled = re.compile(
        r"#\s?<AUTOUPDATE>(?:.*\n)+#\s?<\/AUTOUPDATE>", flags=re.MULTILINE
    )
    if not compiled.search(source):
        err = f"Could not find an existing <AUTOUPDATE> block in {mod.__name__}"
        raise ValueError(err)
    all_block = (
        ALL_TEMPLATE.strip().format(all=",".join(f'"{item}"' for item in all_items))
        if all_items is not None
        else ""
    )
    replaced = compiled.sub(
        SOURCE_TEMPLATE.strip().format(
            all_block=all_block, statements="\n".join(statements or [])
        ),
        source,
    )
    Path(inspect.getfile(mod)).write_text(replaced)


def spec_stub(versions: Iterable[BidsPathSpecFile], latest: BidsPathSpecFile):
    """Generate spec function stub declarations.

    Parameters
    ----------
    versions
        Versions for which to generate declarations
    latest
        Version to which the `latest` stub should point
    """
    for version in versions:
        yield spec_func.format_pyi(version)

    yield spec_func.format_pyi(latest)


def get_latest(versions: Iterable[BidsPathSpecFile]) -> tuple[str, BidsPathSpecFile]:
    """Get the latest version according to semver order."""
    all_versions: dict[str, BidsPathSpecFile] = {}

    for version in versions:
        all_versions[version["version"]] = version

    version = sorted(all_versions, key=lambda v: tuple(v.split(".")))[-1]
    return version.replace(".", "_"), {
        **all_versions[version],
        "version": "latest",
    }


def main():
    """update_bids entrypoint."""
    all_specs = sorted(get_specs(), key=lambda o: o["version"])
    latest_version, latest_spec = get_latest(all_specs)
    generate_stub(
        specs,
        ["from ._utils import BidsPathSpec", "LATEST: str"],
        spec_stub(all_specs, latest_spec),
    )

    versions = [spec["version"].replace(".", "_") for spec in all_specs]
    spec_list = ",".join(f'"{v}"' for v in versions)
    valid_specs = spec_list + ', "latest"'
    update_source(_config, None, [f"VALID_SPECS: TypeAlias = Literal[{valid_specs}]"])

    update_source(
        specs,
        it.chain((version for version in versions), ("latest", "LATEST")),
        [f"_SPECS = [{spec_list}]", f'# LATEST = "{latest_version}"'],
    )


if __name__ == "__main__":
    main()
