""" Minor modification to allow glob_wildcards() to have multiple occurence of
same wildcard """

import collections
import os
import re
from itertools import chain


def regex(filepattern):
    """Build Snakebids regex based on the given file pattern."""
    regex_list = []
    last = 0
    wildcards = set()
    for match in _wildcard_regex.finditer(filepattern):
        regex_list.append(re.escape(filepattern[last : match.start()]))
        wildcard = match.group("name")
        if wildcard in wildcards:
            if match.group("constraint"):
                raise ValueError(
                    "Constraint regex must be defined only in the first "
                    "occurence of the wildcard in a string."
                )
            regex_list.append(f"(?P={wildcard})")
        else:
            wildcards.add(wildcard)
            constraint = (
                match.group("constraint") if match.group("constraint") else ".+"
            )
            regex_list.append(f"(?P<{wildcard}>{constraint})")

        last = match.end()
    regex_list.append(re.escape(filepattern[last:]))
    regex_list.append("$")  # ensure that the match spans the whole file
    return "".join(regex_list)


_wildcard_regex = re.compile(
    r"""
    \{
        (?=(   # This lookahead assertion emulates an 'atomic group'
               # which is required for performance
            \s*(?P<name>\w+)                    # wildcard name
            (\s*,\s*
                (?P<constraint>                 # an optional constraint
                    ([^{}]+ | \{\d+(,\d+)?\})*  # allow curly braces to nest
                )                               # ...  as in '{w,a{3,5}}'
            )?\s*
        ))\1
    \}
    """,
    re.VERBOSE,
)


def glob_wildcards(pattern, files=None, followlinks=False):
    """Glob the values of wildcards by matching a pattern to the filesystem.

    Parameters
    ----------
    pattern : str
        Path including wildcards to glob on the filesystem.
    files : list of str, optional
        Files from which to glob wildcards. If None (default), the directory
        corresponding to the first wildcard in the pattern is walked, and
        wildcards are globbed from all files.
    followlinks : bool, optional
        Whether to follow links when globbing wildcards. Default: False.

    Returns
    -------
    namedtuple
        "Wildcards" named tuple where each field name is the name of a
        wildcard and the value of each field is a list of values for the
        corresponding wildcard.
    """
    pattern = os.path.normpath(pattern)
    first_wildcard = re.search("{[^{]", pattern)
    dirname = (
        os.path.dirname(pattern[: first_wildcard.start()])
        if first_wildcard
        else os.path.dirname(pattern)
    )
    if not dirname:
        dirname = "."

    names = [match.group("name") for match in _wildcard_regex.finditer(pattern)]

    # remove duplicates while preserving ordering
    names = list(collections.OrderedDict.fromkeys(names))

    # TODO: using namedtuple prevents python keywords (as, from, etc) from being used
    #       as wildcards. A Dict would be more appropriate here
    # pylint: disable-msg=invalid-name
    Wildcards = collections.namedtuple("Wildcards", names)

    wildcards = Wildcards(*[[] for _ in names])

    pattern = re.compile(regex(pattern))

    if files is None:
        files = (
            os.path.normpath(os.path.join(dirpath, f))
            for dirpath, dirnames, filenames in os.walk(
                dirname, followlinks=followlinks
            )
            for f in chain(filenames, dirnames)
        )

    for f in files:
        match = re.match(pattern, f)
        if match:
            for name, value in match.groupdict().items():
                getattr(wildcards, name).append(value)
    return wildcards


def update_wildcard_constraints(
    pattern, wildcard_constraints, global_wildcard_constraints
):
    """Update wildcard constraints.

    Parameters
    ----------
    pattern : str
        Pattern on which to update constraints.
    wildcard_constraints : dict
        Dictionary of wildcard:constraint key-value pairs.
    global_wildcard_constraints : dict
        Dictionary of wildcard:constraint key-value pairs.
    """

    def replace_constraint(match):
        name = match.group("name")
        constraint = match.group("constraint")
        newconstraint = wildcard_constraints.get(
            name, global_wildcard_constraints.get(name)
        )
        if name in examined_names:
            return match.group(0)
        examined_names.add(name)
        # Don't override if constraint already set
        if constraint is not None:
            return match.group(0)
        # Only update if a new constraint has actually been set
        if newconstraint is not None:
            return f"{{{name},{newconstraint}}}"
        return match.group(0)

    examined_names = set()
    updated = _wildcard_regex.sub(replace_constraint, pattern)

    return updated
