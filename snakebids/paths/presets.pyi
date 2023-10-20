# This stub file is automatically generated
# It can be updated using::
#
#      poetry run poe update_bids

from pathlib import Path

def bids_v0_0_0(
    root: str | Path | None = None,
    *,
    datatype: str | None = None,
    prefix: str | None = None,
    suffix: str | None = None,
    extension: str | None = None,
    **entities: str | bool,
) -> str:
    """Generate bids or bids-like paths

    Path is compiled based on the 'v0.0.0' spec, with the general form::

        [root]/[sub-{subject}]/[ses-{session}]/
        [prefix]_[sub-{subject}]_[ses-{session}]_[{key}-{val}_ ... ]_[suffix]

    If no arguments are specified, an empty string will be returned.


    Parameters
    ----------
    root
        Root folder to include in the path (e.g. ``results``)
    datatype
        Folder to include after sub-/ses- (e.g. ``anat``, ``dwi`` )
    prefix
        String to prepend to the file name. Useful for injecting custom entities at
        the front of the filename, e.g. ``tpl-{tpl}``
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
    ...

def bids(
    root: str | Path | None = None,
    *,
    datatype: str | None = None,
    prefix: str | None = None,
    suffix: str | None = None,
    extension: str | None = None,
    **entities: str | bool,
) -> str:
    """Generate bids or bids-like paths

    Path is compiled based on the 'latest' spec (currently pointing to 'v0_0_0'), with
    the general form::

        [root]/[sub-{subject}]/[ses-{session}]/
        [prefix]_[sub-{subject}]_[ses-{session}]_[{key}-{val}_ ... ]_[suffix]

    If no arguments are specified, an empty string will be returned.


    Parameters
    ----------
    root
        Root folder to include in the path (e.g. ``results``)
    datatype
        Folder to include after sub-/ses- (e.g. ``anat``, ``dwi`` )
    prefix
        String to prepend to the file name. Useful for injecting custom entities at
        the front of the filename, e.g. ``tpl-{tpl}``
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
    ...
