"""Internal utilities for snakemake template."""

from __future__ import annotations

import string
from collections.abc import Iterator, Mapping, Sequence
from typing import TYPE_CHECKING, Any, Final, overload

import attrs
from typing_extensions import LiteralString, override


@attrs.define(frozen=True)
class _Wildcard:
    """Represents a single wildcard with label, constraint, and formatted output."""

    label: str
    constraint: str

    @property
    def wildcard(self) -> str:
        """Return the brace-wrapped wildcard string."""
        return f"{{{self.label},{self.constraint}}}"

    def __str__(self) -> str:
        return self.wildcard


class SnakemakeWildcards:
    """Define and handle optional wildcard syntax for Snakemake.

    This class encapsulates the rules for defining variable wildcards,
    dummy wildcards, and directory wildcards, as well as special wildcards.

    Parameters
    ----------
    tag : str
        The entity tag name (e.g., "sub", "ses", "run").
    """

    # Special wildcard class attributes (static, tag-independent)
    underscore: Final[_Wildcard] = _Wildcard(
        label="___", constraint=r"^|(?<=\/)|(?<=[^\/])_(?=[^\.])"
    )
    slash: Final[_Wildcard] = _Wildcard(
        label="__d__", constraint=r"^|(?<=\/)|(?<=[^\/])\/"
    )
    datatype: Final[_Wildcard] = _Wildcard(
        label="datatype", constraint=r"(?:(?:^|(?<=\/))[^_\/\-\.]+(?=\/))?"
    )
    suffix: Final[_Wildcard] = _Wildcard(
        label="suffix", constraint=r"(?:(?:^|(?<=\/|_))[^_\/\-\.]+)?"
    )
    extension: Final[_Wildcard] = _Wildcard(
        label="extension", constraint=r"(?:\.[^_\/\-]+$)?"
    )

    # Instance attributes (tag-dependent)
    variable: _Wildcard
    dummy: _Wildcard
    directory: _Wildcard

    def __init__(self, tag: str) -> None:
        """Initialize with an entity tag name.

        Parameters
        ----------
        tag : str
            The entity tag name.
        """
        # subject/session are special: the wildcard label uses the full name, but the
        # BIDS tag prefix in file paths (and therefore in regex constraints) is the
        # short form. This mirrors the same special-casing in SnakemakeFormatter.
        # TODO: a future refactor should unify this mapping across the codebase.
        self._wildcard = {"sub": "subject", "ses": "session"}.get(tag, tag)
        self._tag = {"subject": "sub", "session": "ses"}.get(
            self._wildcard, self._wildcard
        )

        # Initialize instance wildcards with tag-dependent constraints
        self.variable = _Wildcard(
            label=self._wildcard,
            constraint=rf"(?:(?<={self._tag}\-)[^_\/\-\.]+(?=_|\/|\.|$))?",
        )
        self.dummy = _Wildcard(
            label=f"_{self._wildcard}_",
            constraint=rf"(?:(?:^|(?<=\/)|(?<=[^\/])_){self._tag}\-(?=[^_\/\-\.]))?",
        )
        self.directory = _Wildcard(
            label=f"_{self._wildcard}_d_",
            constraint=rf"(?:(?:^|(?<=\/)){self._tag}\-[^_\/\-\.]+\/)?",
        )


class MissingEntityError(KeyError):
    """Raised when special wildcard unaccompanied by associated key."""

    def __init__(self, msg: str, dummy: str, entity: str, *args: object) -> None:
        super().__init__(msg, *args)
        self.dummy = dummy
        self.entity = entity

    @classmethod
    def format_msg(cls, dummy: str, entity: str):
        """Return instance with formatted error message."""
        return cls(
            f"Missing required entity '{entity}' for wildcard '{dummy}'", dummy, entity
        )


class SnakemakeFormatter(string.Formatter):
    """Custom formatter for snakemake templates with optional BIDS entities.

    This formatter handles complex logic of optional BIDS entities, applying
    appropriate defaults for dummy wildcards and managing underscore suppression
    based on context.

    Values for the magic wildcards may be explicitly provided and will be used
    transparently when given.

    To ensure output integrity, no entity wildcard may have ``/`` or ``_`` as the final
    character, except for ``__d__`` and ``___``, which if not set to a blank string may
    ONLY be assigned either ``/`` and ``_`` respectively.

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
    _UNEXPECTED_CLOSE = "unexpected '}' in string"
    _MISSING_CLOSE = "expected '}' before end of string"
    _UNEXPECTED_OPEN = "unexpected '{' in field name"

    @override
    def __init__(self, allow_missing: bool = False) -> None:
        super().__init__()
        self.allow_missing = allow_missing
        self._current_field: str | None = None
        self._underscore: str = ""

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
        """Call base vformat after resetting _underscore."""
        self._underscore = ""
        self._current_field = None
        return super().vformat(format_string, args, kwargs)

    @override
    def parse(  # noqa: PLR0912
        self, format_string: str
    ) -> Iterator[tuple[str, str | None, str | None, str | None]]:
        """Parse format string, stripping constraints and format specifications.

        The function is implemented from scratch in order to avoid the special treatment
        of ``!``, ``:``, and ``[`` by ``string.Formatter.parse()``. It is about 4-5
        times slower than the native implementation, but reasonably well optimized.

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
        next_char = "{"
        j = format_string.find("}")
        j = length if j == -1 else j

        while anchor < length:
            i = format_string.find(next_char, i + 1)

            # point i to the earlier character, j to the further
            if i > j:
                i, j = j, i
            elif i == -1:
                # if end of string, terminate
                if j == length:
                    yield format_string[anchor:], None, None, None
                    return
                i = j
                j = length

            # if final character
            if i + 1 == length:
                if format_string[i] == "{":
                    raise ValueError(self._MISSING_CLOSE)
                raise ValueError(self._UNEXPECTED_CLOSE)

            # if doubled brace
            if format_string[i + 1] == format_string[i]:
                i += 1
                yield format_string[anchor:i], None, None, None
                anchor = i + 1
                next_char = format_string[i]
                self._underscore = "_"
                continue

            # if undoubled } outside of field
            if format_string[i] == "}":
                raise ValueError(self._UNEXPECTED_CLOSE)

            # Otherwise, must be {, so start field
            literal_text = format_string[anchor:i]
            if j == length:
                raise ValueError(self._MISSING_CLOSE)
            if format_string.find("{", i + 1, j) != -1:
                raise ValueError(self._UNEXPECTED_OPEN)

            # take only field name before comma
            comma = format_string.find(",", i + 1, j)
            if comma != -1:
                field_name = format_string[i + 1 : comma]
                # Store full text (with constraint) for allow_missing support
                self._current_field = format_string[i + 1 : j]
            else:
                field_name = format_string[i + 1 : j]
                # Clear any stale constraint from a previous occurrence
                self._current_field = None

            # setup for the next field
            i = j
            j = format_string.find("}", j + 1)
            j = length if j == -1 else j
            anchor = i + 1
            next_char = "{"
            if literal_text:
                self._underscore = (
                    "" if literal_text[-1] in self.UNDERSCORE_SQUELCHERS else "_"
                )
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
        # This wrapper sets the next leading underscore based on the returned value
        # and handles missing keys
        try:
            value = self._get_value(key, args, kwargs)
        except KeyError:
            if self.allow_missing:
                # Return original brace-wrapped wildcard including constraint
                if TYPE_CHECKING:
                    assert isinstance(key, str)
                if key == "__d__" or self._is_directory_wildcard(key):
                    self._underscore = ""
                elif not self._underscore:
                    self._underscore = SnakemakeWildcards.underscore.wildcard
                # if _current_field is "", key must necessarily be ""
                original = self._current_field or key
                return f"{{{original}}}"
            raise

        if value:
            self._underscore = "" if value[-1] in self.UNDERSCORE_SQUELCHERS else "_"
        return value

    def _get_value(
        self, key: str | int, args: Sequence[Any], kwargs: Mapping[str, Any]
    ) -> Any:
        # Rule 0. If it's an integer key, use parent behavior
        if isinstance(key, int):
            result = str(super().get_value(key, args, kwargs))
            self._validate(key, result)
            return result

        # Rule 1. If key is found in kwargs, return it directly
        if key in kwargs:
            result = str(kwargs[key])
            self._validate(key, result)
            return result

        leading_underscore = self._underscore
        trailing_slash = ""

        # Rule 4. Handle ___ special wildcard
        if key == "___":
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
            raise MissingEntityError.format_msg(key, entity) from err

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

        return leading_underscore + result + trailing_slash

    def _validate(self, key: str | int, value: str):
        if key == "__d__":
            if value and value != "/":
                msg = f"__d__ assigned invalid value {value!r}, must be '/' or ''"
                raise ValueError(msg)
            return
        if key == "___":
            if value and value != "_":
                msg = f"___ assigned invalid value {value!r}, must be '_' or ''"
                raise ValueError(msg)
            return

        if isinstance(key, str) and self._is_directory_wildcard(key):
            if value.endswith("/"):
                return
            msg = f"entity {key!r} assigned value {value!r} that does not end in '/'"
            raise ValueError(msg)

        if value.endswith(("_", "/")):
            msg = (
                f"entity {key!r} assigned value {value!r} ending with disallowed "
                f"character {value[-1]!r}. This can result in inconsistent behaviour."
            )
            raise ValueError(msg)

    def _is_directory_wildcard(self, key: str):
        min_dir_length = 3
        return key.startswith("_") and key.endswith("_d_") and len(key) > min_dir_length
