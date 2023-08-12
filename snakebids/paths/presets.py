"""Utilities for converting Snakemake apps to BIDS apps."""
from __future__ import annotations

import functools as ft
import itertools as it
import os
import sys
from pathlib import Path

import more_itertools as itx

from snakebids.paths.specs import v0_0_0


@ft.lru_cache
def _parse_spec(include_subject_dir: bool, include_session_dir: bool):
    spec = v0_0_0(subject_dir=include_subject_dir, session_dir=include_session_dir)
    order: list[str] = []
    dirs: set[str] = set()
    aliases: dict[str, str] = {}

    for entry in spec:
        tag = entry.get("tag", entry["entity"])
        order.append(tag)
        aliases[entry["entity"]] = tag
        if entry.get("dir"):
            dirs.add(tag)

    def parse_entities(entities: dict[str, str | bool]) -> dict[str, str]:
        result: dict[str, str] = {}
        for entity, val in entities.items():
            # strip underscores from keys (needed so that users can use reserved
            # keywords by appending a _)
            stripped = entity.rstrip("_")
            unaliased = aliases.get(stripped, stripped)
            if unaliased in result:
                aliased = itx.nth(aliases, list(aliases.values()).index(unaliased))
                err = (
                    "Long and short names of an entity cannot be used in the same "
                    f"call to bids(): got '{aliased}' and '{unaliased}'"
                )
                raise ValueError(err)
            result[unaliased] = str(val)
        return result

    return order, dirs, parse_entities


def bids(
    root: str | Path | None = None,
    datatype: str | None = None,
    prefix: str | None = None,
    suffix: str | None = None,
    extension: str | None = None,
    **entities: str | bool,
) -> str:
    """Helper function for generating bids paths for snakemake workflows.

    File path is of the form::

        [root]/[sub-{subject}]/[ses-{session]/
        [prefix]_[sub-{subject}]_[ses-{session}]_[{key}-{val}_ ... ]_[suffix]

    Parameters
    ----------
    root : str or Path, default=None
        root folder to include in the path (e.g. 'results')
    datatype : str, default=None
        folder to include after sub-/ses- (e.g. anat, dwi )
    prefix : str, default=None
        string to prepend to the file name (typically not defined, unless you
        want tpl-{tpl}, or a datatype)
    suffix : str, default=None
        bids suffix including extension (e.g. 'T1w.nii.gz')
    subject : str, default=None
        subject to use, for folder and filename
    session : str, default=None
        session to use, for folder and filename
    include_subject_dir : bool, default=True
        whether to include the sub-{subject} folder if subject defined
        (default: True)
    include_session_dir : bool, default=True
        whether to include the ses-{session} folder if session defined
        (default: True)
    **entities : dict, optional
        dictionary of bids entities (e.g. space=T1w for space-T1w)

    Returns
    -------
    str
        bids-like file path

    Examples
    --------

    Below is a rule using bids naming for input and output::

        rule proc_img:
            input: 'sub-{subject}_T1w.nii.gz'
            output: 'sub-{subject}_space-snsx32_desc-preproc_T1w.nii.gz'

    With bids() you can instead use::

         rule proc_img:
            input: bids(subject='{subject}',suffix='T1w.nii.gz')
            output: bids(
                subject='{subject}',
                space='snsx32',
                desc='preproc',
                suffix='T1w.nii.gz'
            )

    Note that here we are not actually using "functions as inputs" in snakemake, which
    would require a function definition with wildcards as the argument, and restrict to
    input/params, but bids() is being used simply to return a string.

    Also note that space, desc and suffix are NOT wildcards here, only {subject} is.
    This makes it easy to combine wildcards and non-wildcards with bids-like naming.

    However, you can still use bids() in a lambda function. This is especially useful if
    your wildcards are named the same as bids entities (e.g. {subject}, {session},
    {task} etc..)::

        rule proc_img:
            input: lambda wildcards: bids(**wildcards,suffix='T1w.nii.gz')
            output: bids(
                subject='{subject}',
                space='snsx32',
                desc='preproc',
                suffix='T1w.nii.gz'
            )

    Or another example where you may have many bids-like wildcards used in your
    workflow::

        rule denoise_func:
            input: lambda wildcards: bids(**wildcards, suffix='bold.nii.gz')
            output: bids(
                subject='{subject}',
                session='{session}',
                task='{task}',
                acq='{acq}',
                desc='denoise',
                suffix='bold.nii.gz'
            )

    In this example, all the wildcards will be determined from the output and passed on
    to bids() for inputs. The output filename will have a 'desc-denoise' flag added to
    it.

    Also note that even if you supply entities in a different order, the entities will
    be ordered based on the OrderedDict defined here. If entities not known are
    provided, they will be just be placed at the end (before the suffix), in the order
    you provide them in.

    Notes
    -----

    * For maximum flexibility all arguments are optional (if none are specified, will
      return empty string). Note that datatype and prefix may not be used in isolation,
      but must be given with another entity.

    * Some code adapted from mne-bids, specifically
      https://mne.tools/mne-bids/stable/_modules/mne_bids/utils.html
    """
    if not any([entities, suffix, extension]) and any([datatype, prefix]):
        raise ValueError(
            "At least one of suffix, extension, or an entity must be "
            "supplied.\n\tGot only: "
            + " and ".join(
                filter(
                    None,
                    (
                        f"datatype='{datatype}'" if datatype else None,
                        f"prefix='{prefix}'" if prefix else None,
                    ),
                )
            )
        )

    include_subject_dir = bool(entities.pop("include_subject_dir", True))
    include_session_dir = bool(entities.pop("include_session_dir", True))

    order, dirs, parse_entities = _parse_spec(
        include_subject_dir=include_subject_dir, include_session_dir=include_session_dir
    )
    parsed = parse_entities(entities)

    spec_parts: list[str] = []
    custom_parts: list[str] = []
    split: int = sys.maxsize + 1
    path_parts: list[str] = []

    if root:
        path_parts.append(str(root))
    if prefix:
        spec_parts.append(prefix)
    for entity in order:
        # Check for `*` first so that if user specifies an entity called `*` we don't
        # skip setting the split
        if entity == "*":
            split = len(path_parts)
        elif value := parsed.pop(entity, None):
            spec_parts.append(f"{entity}-{value}")
            if entity in dirs:
                path_parts.append(f"{entity}-{value}")
    for key, value in parsed.items():
        custom_parts.append(f"{key}-{value}")

    if datatype:
        path_parts.append(datatype)
    path_parts.append(
        "_".join(it.chain(spec_parts[:split], custom_parts, spec_parts[split:]))
    )

    tail = f"_{suffix}{extension or ''}" if suffix else extension or ""

    return os.path.join(*path_parts) + tail
