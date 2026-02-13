"""Internal utilities for Snakemake wildcard handling."""

from __future__ import annotations

from typing import Final

from snakebids.utils.utils import BidsEntity


class SnakemakeWildcards:
    """Define and handle optional wildcard syntax for Snakemake.

    This class encapsulates the rules for defining variable wildcards,
    dummy wildcards, and directory wildcards, as well as special wildcards.

    Parameters
    ----------
    tag : str
        The entity tag name (e.g., "sub", "ses", "run").
    """

    # Special wildcard class attributes with their constraints
    # Format: NAME,CONSTRAINT
    underscore: Final[str] = r"___,^|(?<=\/)|(?<!\/)_(?=[^\.])"
    d: Final[str] = r"__d__,^|(?<=\/)|(?<=.)\/"
    datatype: Final[str] = r"datatype,(?:(?:^|(?<=\/))[^_\/\-\n]+(?=\/))?"
    suffix: Final[str] = r"suffix,(?:(?:^|(?<=\/|_))[^_\/\-]+)?"
    extension: Final[str] = r"extension,(?:\.[^_\/\-]+$)?"

    def __init__(self, tag: str) -> None:
        """Initialize with an entity tag name.

        Parameters
        ----------
        tag : str
            The entity tag name.
        """
        # Use BidsEntity to handle tag to wildcard conversion
        self._tag = BidsEntity.from_tag(tag).wildcard

    @property
    def dummy(self) -> str:
        """Return the dummy wildcard in NAME,CONSTRAINT format.

        Returns
        -------
        str
            Dummy wildcard format: _TAG_,CONSTRAINT
        """
        constraint = rf"(?:(?:^|(?<=\/)|(?<!\/)_){self._tag}\-(?=[^_\/\-\n]))?"
        return f"_{self._tag}_,{constraint}"

    @property
    def variable(self) -> str:
        """Return the variable wildcard in NAME,CONSTRAINT format.

        Returns
        -------
        str
            Variable wildcard format: TAG,CONSTRAINT
        """
        constraint = rf"(?:(?<={self._tag}\-)[^_\/\-\n]+(?=_|\/|$))?"
        return f"{self._tag},{constraint}"

    @property
    def directory(self) -> str:
        """Return the directory wildcard in NAME,CONSTRAINT format.

        Returns
        -------
        str
            Directory wildcard format: _TAG_d_,CONSTRAINT
        """
        constraint = rf"(?:(?:^|(?<=\/)){self._tag}\-[^_\/\-\n]+\/)?"
        return f"_{self._tag}_d_,{constraint}"
