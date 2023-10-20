import textwrap

import docstring_parser as docstr

TEMPLATE = '''
def bids{spec}(
    root: str | Path | None = None,
    *,
    datatype: str | None = None,
    prefix: str | None = None,
    suffix: str | None = None,
    extension: str | None = None,
    **entities: str | bool,
) -> str:
    """{docstring}
    """
    ...
'''
DOCSTRING = """Generate bids or bids-like paths

Path is compiled based on the '{spec_label}' spec{spec_clarify}, with the
general form::

    [root]/[sub-{{subject}}]/[ses-{{session}}]/
    [prefix]_[sub-{{subject}}]_[ses-{{session}}]_[{{key}}-{{val}}_ ... ]_[suffix]

If no arguments are specified, an empty string will be returned.


Parameters
----------

root
    Root folder to include in the path (e.g. ``results``)
datatype
    Folder to include after sub-/ses- (e.g. ``anat``, ``dwi`` )
prefix
    String to prepend to the file name. Useful for injecting custom entities at
    the front of the filename, e.g. ``tpl-{{tpl}}``
suffix
    Suffix plus, optionally, the extension (e.g. ``T1w.nii.gz``)
extension
    bids extension, beginning with ``.`` (e.g. ``.nii.gz``).  Typically
    shouldn't be specified manually: extensions should be listed along with the
    suffix.
entities
    bids entities as keyword arguments paired with values (e.g. ``space="T1w"``
    for ``space-T1w``)

"""


def format_pyi(spec: str, spec_label: str, spec_clarify: str = ""):
    doc = docstr.parse(
        DOCSTRING.format(spec_label=spec_label, spec_clarify=spec_clarify)
    )
    if doc.long_description:
        doc.long_description = "\n\n".join(
            "\n".join(textwrap.wrap(para, 84)) if not para.startswith("    ") else para
            for para in doc.long_description.split("\n\n")
        )
    return TEMPLATE.format(
        spec=spec,
        docstring=docstr.compose(doc),
    )
