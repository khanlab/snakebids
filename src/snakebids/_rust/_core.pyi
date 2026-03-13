"""Type stubs for the optional compiled Rust extension ``snakebids._rust._core``."""

def parse_format_string(
    format_string: str,
) -> list[tuple[str, str | None, bool | None, str | None]]:
    """Parse a Snakemake-style format string character by character.

    Parameters
    ----------
    format_string : str
        The format string to parse.

    Returns
    -------
    list[tuple[str, str | None, bool | None, str | None]]
        A list of ``(literal_text, field_name, squelch_underscore, constraint)`` tuples:

        - ``literal_text``       – literal text before the next field
        - ``field_name``         – ``None`` for non-field segments; the field name otherwise
        - ``squelch_underscore`` – ``None`` = no change to ``_underscore`` (empty literal or
          end-of-string); ``True`` = set ``_underscore = ""``; ``False`` = set ``_underscore = "_"``
        - ``constraint``         – ``None`` for non-field segments; ``""`` for a field with no
          constraint; ``",<pattern>"`` for a field with a constraint; concatenating
          ``field_name + constraint`` gives the full ``_current_field`` value

    Raises
    ------
    ValueError
        On malformed input (unexpected ``}`` / ``{`` or missing ``}``).
    """
    ...
