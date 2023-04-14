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
