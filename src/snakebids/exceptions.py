from __future__ import annotations

from collections.abc import Iterable


class ConfigError(Exception):
    """Exception raised for errors with the Snakebids config."""

    def __init__(self, msg: str) -> None:
        self.msg = msg
        super().__init__(msg)


class RunError(Exception):
    """Exception raised for errors in generating and running the snakemake workflow."""

    def __init__(self, msg: str, *args: object) -> None:
        super().__init__(msg, *args)
        self.msg = msg


class PybidsError(Exception):
    """Exception raised when pybids encounters a problem."""


class DuplicateComponentError(Exception):
    """Raised when a dataset is constructed from components with the same name."""

    def __init__(self, duplicated_names: Iterable[str]):
        self.duplicated_names_str = ", ".join(map(repr, duplicated_names))
        super().__init__(
            "A BidsDataset cannot be instantiated from BidsComponents with the same "
            "names. The following duplicate names were found: "
            f"{self.duplicated_names_str}."
        )


class SnakebidsPluginError(Exception):
    """Exception raised when a Snakebids plugin encounters a problem."""
