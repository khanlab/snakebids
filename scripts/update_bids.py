"""Recompile the bids function stub file based on latest specs"""
from __future__ import annotations

import inspect
import itertools as it
import re
from pathlib import Path
from types import ModuleType
from typing import Iterable

import black

from snakebids.paths import presets, specs
from snakebids.paths._templates import bids_func, spec_func
from snakebids.paths.utils import BidsPathSpecFile, get_specs


def generate_stub(mod: ModuleType, imports: list[str], funcs: Iterable[str]):
    prelude = [
        "# This stub file is automatically generated",
        "# It can be updated using::",
        "#",
        "#      poetry run poe update_bids",
        "",
    ]
    assert mod.__file__
    # pyi = "\n".join(it.chain(prelude, imports, funcs))
    # print(pyi)
    # print(black.format_str(
    #     pyi,
    #     mode=black.Mode(is_pyi=True),
    # ))
    with Path(mod.__file__).with_suffix(".pyi").open("w") as f_:
        f_.write(
            black.format_str(
                "\n".join(it.chain(prelude, imports, funcs)),
                mode=black.Mode(is_pyi=True),
            )
        )


SOURCE_TEMPLATE = """
# <AUTOUPDATE>
# The code between these tags is automatically generated. Do not
# manually edit
# To update, run::
#
#       poetry run poe update_bids
#

if not TYPE_CHECKING:
    __all__ = [ # noqa:F822
        {all}
    ]

    def __dir__():
        return __all__

{statements}
# </AUTOUPDATE>
"""


def update_source(
    mod: ModuleType, all_items: Iterable[str], statements: list[str] | None = None
):
    source = inspect.getsource(mod)
    all_formatted = ",".join(f'"{item}"' for item in all_items)
    compiled = re.compile(
        r"#\s?<AUTOUPDATE>(?:.*\n)+#\s?<\/AUTOUPDATE>", flags=re.MULTILINE
    )
    if not compiled.search(source):
        err = f"Could not find an existing <AUTOUPDATE> block in " f"{mod.__name__}"
        raise ValueError(err)
    replaced = black.format_str(
        compiled.sub(
            SOURCE_TEMPLATE.strip().format(
                all=all_formatted, statements="\n".join(statements or [])
            ),
            source,
        ),
        mode=black.Mode(),
    )
    with open(inspect.getfile(mod), "w") as f_:
        f_.write(replaced)


def presets_stub(versions: Iterable[str], latest: str):
    for member in versions:
        yield bids_func.format_pyi(
            spec=f"_{member.replace('.', '_')}", spec_label=member
        )

    yield bids_func.format_pyi(
        spec="",
        spec_label="latest",
        spec_clarify=f" (currently pointing to '{latest}')",
    )


def spec_stub(versions: Iterable[BidsPathSpecFile], latest: BidsPathSpecFile):
    for version in versions:
        yield spec_func.format_pyi(version)

    yield spec_func.format_pyi(latest)


def get_latest(versions: Iterable[BidsPathSpecFile]) -> tuple[str, BidsPathSpecFile]:
    all_versions: dict[str, BidsPathSpecFile] = {}

    for version in versions:
        all_versions[version["version"]] = version

    version = sorted(all_versions, key=lambda v: tuple(v.split(".")))[-1]
    return version.replace(".", "_"), {
        **all_versions[version],
        "version": "latest",
    }


def main():
    all_specs = list(get_specs())
    latest_version, latest_spec = get_latest(all_specs)
    generate_stub(
        specs,
        ["from .utils import BidsPathSpec", "LATEST: str"],
        spec_stub(all_specs, latest_spec),
    )
    generate_stub(
        presets,
        ["from pathlib import Path"],
        presets_stub((spec["version"] for spec in all_specs), latest_version),
    )

    versions = [spec["version"].replace(".", "_") for spec in all_specs]
    update_source(
        presets, it.chain((f"bids_{version}" for version in versions), ("bids",))
    )

    spec_list = ",".join(f'"{v}"' for v in versions)
    update_source(
        specs,
        it.chain((version for version in versions), ("latest", "LATEST")),
        [f"_SPECS = [{spec_list}]", f'LATEST = "{latest_version}"'],
    )


if __name__ == "__main__":
    main()
