"""Internal utilities for snakemake template."""

from __future__ import annotations

import string
from collections.abc import Iterator, Mapping, Sequence
from typing import Any, Final, overload

from typing_extensions import LiteralString, override

from snakebids.utils.utils import BidsEntity


class SnakemakeFormatter(string.Formatter):
    """Custom formatter for snakemake templates with optional BIDS entities.

    This formatter handles complex logic of optional BIDS entities, applying
    appropriate defaults for dummy wildcards and managing underscore suppression
    based on context.

    Wildcard Label Syntax:
    - Ordinary wildcards: {TAG} (e.g., {subject}, {acq}, {run})
    - Dummy wildcards: {_TAG_} (e.g., {_subject_}, {_acq_}, {_run_})
    - Directory wildcards: {_TAG_d_} (e.g., {_subject_d_}, {_session_d_})
    - Special wildcards: {__d__}, {___}, {datatype}, {suffix}, {extension}

    When formatting, entities should generally be specified using their tag names, e.g.
    ``acq`` instead of ``acquisition``. The two exceptions are ``subject`` and
    ``session``, which are recognized and automatically converted.
    """

    UNDERSCORE_SQUELCHERS: Final = {"/", "_"}

    @override
    def __init__(self) -> None:
        super().__init__()
        self.squelch_underscore = True

    @overload
    def vformat(
        self,
        format_string: LiteralString,
        args: Sequence[LiteralString],
        kwargs: Mapping[LiteralString, LiteralString],
    ) -> LiteralString: ...
    @overload
    def vformat(
        self, format_string: str, args: Sequence[Any], kwargs: Mapping[Any, Any]
    ) -> str: ...

    def vformat(
        self, format_string: str, args: Sequence[Any], kwargs: Mapping[Any, Any]
    ) -> str:
        """Call base vformat after resetting squelch_underscore."""
        self.squelch_underscore = True
        return super().vformat(format_string, args, kwargs)

    @override
    def parse(
        self, format_string: str
    ) -> Iterator[tuple[str, str | None, str | None, str | None]]:
        """Parse format string, stripping constraints and format specifications.

        The function is implemented from scratch in order to avoid the special treatment
        of ``!``, ``:``, and ``[`` by ``string.Formatter.parse()``. It is about 4-5
        times slower than the native implementation, but reasonably well optimized.

        Unlike the native implementation, no error is raised when a singleton ``}`` is
        encountered outside of a field. This is primarily for efficiency, but also
        improves compatibility with snakemake templates, which allow such ``}``.

        Specifying conversion and format strings is not possible. Any text following
        ``,`` is treated as a constraint and discarded. Snakemake constraints allow
        a single level of nested brace pairs ``{}``. This is not supported by this
        formatter.

        Parameters
        ----------
        format_string : str
            The format string to parse

        Yields
        ------
            Each iteration yields (literal_text, field_name, ""|None, None)
        """
        i = -1
        anchor = 0
        length = len(format_string)
        while anchor < length:
            i = format_string.find("{", i + 1)
            if i < 0:
                yield format_string[anchor:], None, None, None
                return
            if i + 1 == length:
                msg = "expected '}' before end of string"
                raise ValueError(msg)
            if format_string[i + 1] == "{":
                i += 1
                self.squelch_underscore = False
                yield format_string[anchor:i], None, None, None
                anchor = i + 1
                continue
            literal_text = format_string[anchor:i]
            field_name = None
            close = format_string.find("}", i + 1)
            if close < 0:
                msg = "expected '}' before end of string"
                raise ValueError(msg)
            if format_string.find("{", i + 1, close) != -1:
                msg = "unexpected '{' in field name"
                raise ValueError(msg)
            if (comma := format_string.find(",", i + 1, close)) != -1:
                field_name = format_string[i + 1 : comma]
            else:
                field_name = format_string[i + 1 : close]

            i = close
            anchor = i + 1
            if literal_text:
                self.squelch_underscore = literal_text[-1] in self.UNDERSCORE_SQUELCHERS
            yield literal_text, field_name, "", None

    @override
    def get_value(
        self, key: str | int, args: Sequence[Any], kwargs: Mapping[str, Any]
    ) -> Any:
        """Get value for wildcard using serial parsing strategy.

        The following rules apply:

        1. If key is found in kwargs, return it directly.
        2. If key is a directory (_TAG_d_) or dummy (_TAG_) type label, return default
           if TAG given in kwargs, or blank string if TAG given as blank string,
           otherwise error.
        3. If ``key == "__d__"``, treat as above, using ``datatype`` as TAG.
        4. If ``key == "___"`` return ``_`` if it does not immediately follow ``/`` or
           ``_``.

        The default for directory type labels is ``TAG-VALUE/``. For dummy type, it is
        ``_TAG-``, excluding the initial ``_`` if the wildcard is at the beginning of
        the string or immediately follows ``/`` or ``_``. For ``__d__``, the default
        is ``/``.

        More generally, the idea is to prevent two adjacent underscores, an underscore
        following a slash, or an underscore at the beginning of the string through
        the filling of structural wildcards.

        Parameters
        ----------
        key : str | int
            The field name to retrieve
        args : tuple[Any, ...]
            Positional arguments (unused)
        kwargs : dict[str, Any]
            Keyword arguments containing entity values

        Returns
        -------
        Any
            The formatted value for the wildcard

        Raises
        ------
        KeyError
            When a required entity is missing from kwargs
        """
        # Rule 0. If it's an integer key, use parent behavior
        if isinstance(key, int):
            return super().get_value(key, args, kwargs)

        # Rule 1. If key is found in kwargs, return it directly
        if key in kwargs:
            # No support for conversion bits, so formatting will be as string
            value = str(kwargs[key])
            # Update squelch_underscore based on returned value
            if value:
                self.squelch_underscore = value[-1] in self.UNDERSCORE_SQUELCHERS
            return value

        leading_underscore = "" if self.squelch_underscore else "_"
        trailing_slash = ""

        # Rule 4. Handle ___ special wildcard
        if key == "___":
            self.squelch_underscore = True
            return leading_underscore

        # Rule 3. Handle __d__ using datatype as entity
        if key == "__d__":
            entity = "datatype"
            trailing_slash = "/"
            leading_underscore = ""

        # Rule 2. Handle directory and dummy wildcards
        elif key.startswith("_") and key.endswith("_") and len(key) > 1:
            entity = key[1:-1]

        # Otherwise, not a special key, so error
        else:
            raise KeyError(key)

        # Rule 2.5. Handle directory wildcards
        if is_directory := entity.endswith("_d"):
            entity = entity[:-2]
            trailing_slash = "/"
            leading_underscore = ""

        # Get provided value for entity
        try:
            entity_value = str(kwargs[entity])
        except KeyError as err:
            msg = f"Missing required entity '{entity}' for wildcard '{key}'"
            raise KeyError(msg) from err

        # Return nothing if value blank
        if not entity_value:
            return ""

        # Create the tag (datatype has no tag)
        result = {"subject": "sub-", "session": "ses-", "datatype": ""}.get(
            entity, entity + "-"
        )

        # If a directory type, append the value
        if is_directory:
            result += entity_value

        self.squelch_underscore = bool(trailing_slash)

        return leading_underscore + result + trailing_slash


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
    underscore: Final[str] = r"___,^|(?<=\/)|(?<=[^\/])_(?=[^\.])"
    slash: Final[str] = r"__d__,^|(?<=\/)|(?<=.)\/"
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
        constraint = rf"(?:(?:^|(?<=\/)|(?<=[^\/])_){self._tag}\-(?=[^_\/\-\n]))?"
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
