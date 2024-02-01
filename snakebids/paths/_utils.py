from __future__ import annotations

from typing import Iterator, List

import importlib_resources as impr
import more_itertools as itx
from typing_extensions import NotRequired, TypeAlias, TypedDict

from snakebids.io.yaml import get_yaml_io
from snakebids.paths import resources

__all__ = ["BidsPathEntitySpec", "BidsPathSpec", "BidsPathSpecFile"]


class BidsPathEntitySpec(TypedDict):
    """Defines an entity in a bids path."""

    entity: str
    """Entity full name"""

    tag: NotRequired[str]
    """Short entity name, as appears in the path"""

    dir: NotRequired[bool]
    """If true, a directory with the entity-value pair is created"""


BidsPathSpec: TypeAlias = List[BidsPathEntitySpec]
"""List of :class:`BidsPathEntitySpec`, defining the order of entities in a bids path"""


class BidsPathSpecFile(TypedDict):
    """Defines the valid structure for a BidsSpec definition file."""

    version: str
    """Version of the spec in semver"""

    description: NotRequired[str]
    """Optional description to be used in the spec function docstring"""

    spec: BidsPathSpec
    """Definition of the spec"""


def load_spec(path: impr.abc.Traversable) -> BidsPathSpecFile:
    """Return the spec saved in the provided yaml file."""
    return get_yaml_io().load(path.read_text())


def get_specs() -> Iterator[BidsPathSpecFile]:
    """Yield all defined bids spec file objects."""
    for path in impr.files(resources).iterdir():
        if (
            path.is_file()
            and path.name.startswith("spec.")
            and path.name.endswith(".yaml")
        ):
            yield load_spec(path)


def get_spec_path(version: str) -> impr.abc.Traversable:
    """Return path corresponding to provided spec version attribute.

    Parameters
    ----------
    version
        version attribute formatted as vx_x_x where x_x_x is spec semver
    """
    dotted = version[1:].replace("_", ".")
    return impr.files(resources).joinpath(f"spec.{dotted}.yaml")


def find_entity(spec: BidsPathSpec, entity: str) -> BidsPathEntitySpec:
    """Return configuration for specified entity out of BidsPathSpec."""
    return itx.one(item for item in spec if item["entity"] == entity)
