"""Internal utilities for Snakemake wildcard handling."""

from __future__ import annotations


class SnakemakeWildcards:
    """Define and handle optional wildcard syntax for Snakemake.

    This class encapsulates the rules for defining variable wildcards,
    dummy wildcards, and directory wildcards, as well as special wildcards.

    Note
    ----
    This class is an internal implementation detail and should not be mentioned
    in public documentation.

    Parameters
    ----------
    tag : str
        The entity tag name (e.g., "sub", "ses", "run"). The tag "sub" will be
        replaced with "subject" and "ses" with "session".
    """

    # Special wildcard class attributes with their constraints
    # Format: NAME,CONSTRAINT (without braces)
    underscore = "___,^|(?<=\\/)|(?<!\\/)_(?=[^\\.])"
    d = "__d__,^|(?<=\\/)|(?<=.)\\/"
    datatype = "datatype,(?:(?:^|(?<=\\/))[^_\\/\\-\\n]+(?=\\/))?"
    suffix = "suffix,(?:(?:^|(?<=\\/|_))[^_\\/\\-]+)?"
    extension = "extension,(?:\\.[^_\\/\\-]+$)?"

    def __init__(self, tag: str) -> None:
        """Initialize with an entity tag name.

        Parameters
        ----------
        tag : str
            The entity tag name. "sub" will be replaced with "subject" and
            "ses" with "session".
        """
        # Replace sub/ses with subject/session
        if tag == "sub":
            self._tag = "subject"
        elif tag == "ses":
            self._tag = "session"
        else:
            self._tag = tag

    @property
    def dummy(self) -> str:
        """Return the dummy wildcard in NAME,CONSTRAINT format.

        Returns
        -------
        str
            Dummy wildcard format: _TAG_,CONSTRAINT (without braces)
        """
        constraint = f"(?:(?:^|(?<=\\/)|(?<!\\/)_){self._tag}\\-(?=[^_\\/\\-\\n]))?"
        return f"_{self._tag}_,{constraint}"

    @property
    def variable(self) -> str:
        """Return the variable wildcard in NAME,CONSTRAINT format.

        Returns
        -------
        str
            Variable wildcard format: TAG,CONSTRAINT (without braces)
        """
        constraint = f"(?:(?<={self._tag}\\-)[^_\\/\\-\\n]+(?=_|\\/|$))?"
        return f"{self._tag},{constraint}"

    @property
    def directory(self) -> str:
        """Return the directory wildcard in NAME,CONSTRAINT format.

        Returns
        -------
        str
            Directory wildcard format: _TAG_d_,CONSTRAINT (without braces)
        """
        constraint = f"(?:(?:^|(?<=\\/)){self._tag}\\-[^_\\/\\-\\n]+\\/)?"
        return f"_{self._tag}_d_,{constraint}"
