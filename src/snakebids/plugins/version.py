from __future__ import annotations

import importlib.metadata as impm
from typing import Any, cast

import attrs

from snakebids import bidsapp
from snakebids.plugins.base import PluginBase


@attrs.define(kw_only=True)
class Version(PluginBase):
    """Expose app version in config.

    A version string can either be explicitly provided, or the name of the distribution.
    In the latter case, :mod:`importlib.metadata` will be used to look up the version.
    Note that this method does not work when the app is not installed; in this case, the
    version will be set to ``unknown``.

    Parameters
    ----------
    version
        Explicit version string
    distribution
        Name of distribution (as installed with ``pip``) for version lookup

    Raises
    ------
    ValueError
        If neither version nor distribution are provided, or if both are provided
    """

    version: str | None = None
    distribution: str | None = None

    PREFIX = "plugins.version"

    def __attrs_post_init__(self):
        if self.version is not None and self.distribution is not None:
            msg = "version and distribution may not both be specified"
            raise ValueError(msg)
        if self.version is None and self.distribution is None:
            msg = "One of version or distribution must be specified"
            raise ValueError(msg)

    @bidsapp.hookimpl
    def initialize_config(self, config: dict[str, Any]):
        """Assign version to config."""
        if self.version is not None:
            version = self.version
        else:
            try:
                version = impm.version(cast(str, self.distribution))
            except impm.PackageNotFoundError:
                version = "unknown"
        self.assign(config, "version", version)
