from __future__ import annotations

from pathlib import Path

from snakebids.paths._config import get_bids_func
from snakebids.paths.specs import LATEST


def bids(
    root: str | Path | None = None,
    *,
    datatype: str | None = None,
    prefix: str | None = None,
    suffix: str | None = None,
    extension: str | None = None,
    **entities: str | bool,
) -> str:
    """Generate bids or bids-like paths.

    File path is of the form::

        [root]/[sub-{subject}]/[ses-{session]/
        [prefix]_[sub-{subject}]_[ses-{session}]_[{key}-{val}_ ... ]_[suffix]

    If no arguments are specified, an empty string will be returned.

    Datatype and prefix may not be used in isolation, but must be given with
    another entity.

    BIDS paths are built based on specs, which are versioned for long-term stability.
    The latest version is ``<version>``. Information on its spec can be found at
    :func:`~snakebids.paths.specs.<version>`.

    .. warning::

        By default, ``bids()`` will always use the latest BIDS spec. This is unsafe for
        production environments, as the spec may be updated without warning, even on
        patch releases. These updates may change the path output by ``bids()``,
        resulting in breaking changes in downstream apps

        Production code should always explicitly set the spec version using
        :func:`~snakebids.set_bids_spec`::

            from snakebids import set_bids_spec
            set_bids_spec("<version>")

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


    Examples
    --------
    Below is a rule using bids naming for input and output::

        rule proc_img:
            input: 'sub-{subject}_T1w.nii.gz' output:
            'sub-{subject}_space-snsx32_desc-preproc_T1w.nii.gz'

    With bids() you can instead use::

            rule proc_img: input: bids(subject='{subject}',suffix='T1w.nii.gz')
            output: bids(
                subject='{subject}', space='snsx32', desc='preproc',
                suffix='T1w.nii.gz'
            )

    Note that here we are not actually using "functions as inputs" in snakemake,
    which would require a function definition with wildcards as the argument, and
    restrict to input/params, but bids() is being used simply to return a string.

    Also note that space, desc and suffix are NOT wildcards here, only {subject} is.
    This makes it easy to combine wildcards and non-wildcards with bids-like naming.

    However, you can still use bids() in a lambda function. This is especially
    useful if your wildcards are named the same as bids entities (e.g. {subject},
    {session}, {task} etc..)::

        rule proc_img:
            input: lambda wildcards: bids(**wildcards,suffix='T1w.nii.gz') output:
            bids(
                subject='{subject}', space='snsx32', desc='preproc',
                suffix='T1w.nii.gz'
            )

    Or another example where you may have many bids-like wildcards used in your
    workflow::

        rule denoise_func:
            input: lambda wildcards: bids(**wildcards, suffix='bold.nii.gz') output:
            bids(
                subject='{subject}', session='{session}', task='{task}',
                acq='{acq}', desc='denoise', suffix='bold.nii.gz'
            )

    In this example, all the wildcards will be determined from the output and passed
    on to bids() for inputs. The output filename will have a 'desc-denoise' flag
    added to it.

    Also note that even if you supply entities in a different order, the entities
    will be ordered based on the OrderedDict defined here. If entities not known are
    provided, they will be just be placed at the end (before the suffix), in the
    order you provide them in.
    """
    return get_bids_func()(
        root,
        datatype=datatype,
        prefix=prefix,
        suffix=suffix,
        extension=extension,
        **entities,
    )


assert bids.__doc__  # noqa: S101
bids.__doc__ = bids.__doc__.replace("<version>", LATEST)
