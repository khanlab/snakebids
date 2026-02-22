from __future__ import annotations

from collections.abc import Iterator
from enum import Enum, auto
from typing import TypeAlias

import importlib_resources as impr
import more_itertools as itx
from typing_extensions import NotRequired, TypedDict

from snakebids.io.yaml import get_yaml_io
from snakebids.paths import resources

__all__ = [
    "OPTIONAL_WILDCARD",
    "BidsFlags",
    "BidsPathEntitySpec",
    "BidsPathSpec",
    "BidsPathSpecFile",
]


class BidsPathEntitySpec(TypedDict):
    """Defines an entity in a bids path."""

    entity: str
    """Entity full name"""

    tag: NotRequired[str]
    """Short entity name, as appears in the path"""

    dir: NotRequired[bool]
    """If true, a directory with the entity-value pair is created"""


BidsPathSpec: TypeAlias = list[BidsPathEntitySpec]
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


class BidsFlags(Enum):
    """Enum for indicating special template behaviors in Snakemake wildcards.

    This enum is used to specify optional entities in BidsComponent wildcards.
    When used as a value in a wildcards dictionary, it signals that the
    entity should be treated as optional in the generated Snakemake template.
    """

    OPTIONAL_WILDCARD = auto()
    """Indicates that an entity should be treated as optional in Snakemake templates.

    When this value is used for an entity in a wildcards dictionary, the `bids()`
    function will generate a constrained wildcard that allows the entity to be
    present or absent in the path.

    Example
    -------
    >>> from snakebids.paths import bids, OPTIONAL_WILDCARD
    >>> # Generate a path with optional 'acq' entity
    >>> template = bids(
    ...     subject="{subject}",
    ...     session="{session}",
    ...     acq=BidsFlags.OPTIONAL_WILDCARD,
    ...     suffix="T1w.nii.gz"
    ... )

    .. warning::
        Templates with optional wildcards are not compatible with Python's
        `str.format()` method. They use Snakemake-specific constraint syntax
        that can only be interpreted by Snakemake's wildcard resolution system.
    """


OPTIONAL_WILDCARD = BidsFlags.OPTIONAL_WILDCARD
