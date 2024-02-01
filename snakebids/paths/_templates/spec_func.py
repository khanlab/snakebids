from __future__ import annotations

import textwrap
from typing import TYPE_CHECKING

from snakebids.utils.utils import entity_to_wildcard

if TYPE_CHECKING:
    from snakebids.paths._utils import BidsPathSpec, BidsPathSpecFile

TEMPLATE = '''
def {spec}(subject_dir: bool = True, session_dir: bool = True) -> BidsPathSpec:
    """{docstring}
    """
    ...
'''
DOCSTRING = """{description}

Parameters
----------
subject_dir
    If False, downstream path generator will not include the subject dir
    `sub-{{subject}}/*`
session_dir : bool, optional
    If False, downstream path generator will not include the session dir
    `*/ses-{{session}}/*`

"""

DEFAULT_DESCRIPTION = """Bids Spec {version}

Supply this to snakebids.bids_factory to construct a corresponding bids function
"""


def _wrap_template(template: str, length: int):
    def recurse(lines: list[str]) -> list[str]:
        line = lines[-1]
        if len(line) <= length:
            return lines
        i = line[length - 1 :: -1].index("_")
        return recurse(lines[:-1] + [line[: length - i], "    " + line[length - i :]])

    return "\n".join(recurse([template]))


def compile_example(spec: BidsPathSpec):
    # import within function to avoid circular import
    from snakebids.paths._factory import bids_factory

    entities = [listing["entity"] for listing in spec]
    standard_entities = ("prefix", "datatype", "suffix", "extension")
    try:
        wild = entities.index("*")
    except ValueError:
        wild = -2
    template = bids_factory(spec)(
        **entity_to_wildcard(standard_entities),
        **entity_to_wildcard(e for e in entities if e != "*"),
    )
    search = f"{spec[wild+1].get('tag', entities[wild+1])}-{{{entities[wild+1]}}}_"
    i = template.index(search)
    if wild < 0:
        i = len(search) + i
    return _wrap_template(template[:i] + "..._" + template[i:], 80)


def _import_docstring_parser():
    """Isolated import function that can be mocked in tests."""
    import docstring_parser as docstr

    return docstr


def format_doc(spec: BidsPathSpecFile):
    if (description := spec.get("description")) is None:
        description = DEFAULT_DESCRIPTION.format(version=spec["version"])

    try:
        docstr = _import_docstring_parser()
    except ImportError:
        return description.strip()

    doc = docstr.parse(DOCSTRING.format(description=description.strip()))
    if doc.long_description:
        doc.long_description = (
            "\n\n".join(
                "\n".join(textwrap.wrap(para, 84))
                if not para.startswith("    ")
                else para
                for para in doc.long_description.split("\n\n")
            )
            + "\n\nFormatted as::\n\n    "
            + compile_example(spec["spec"])
        )

    return docstr.compose(doc)


def format_pyi(spec: BidsPathSpecFile):
    return TEMPLATE.format(
        spec=spec["version"].replace(".", "_"),
        docstring=format_doc(spec),
    )
