def version(distribution_name: str) -> str:
    """Get the version string for the named package.

    :param distribution_name: The name of the distribution package to query.
    :return: The version string for the package as defined in the package's
        "Version" metadata key.
    """
    ...

class PackageNotFoundError(ModuleNotFoundError):
    """The package was not found."""

    def __str__(self) -> str: ...
    @property
    def name(self) -> str: ...
