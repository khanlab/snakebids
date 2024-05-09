from argparse import ArgumentParser
from pathlib import Path
from types import ModuleType
from typing import Any, Callable, Iterable, Sequence

from snakemake.common import configfile as configfile  # type: ignore

configfile: ModuleType

class WildcardError(Exception): ...

def load_configfile(configpath: str) -> dict[str, Any]:
    "Load a JSON or YAML configfile as a dict, then checks that it's a dict."

def expand(
    filepatterns: Sequence[Path | str] | Path | str,
    func: Callable[[Iterable[str]], Iterable[Iterable[str]]] | None = ...,
    /,
    *,
    allow_missing: bool | Sequence[str] | str = ...,
    **wildcards: Sequence[str] | str,
) -> list[str]:
    """
    Expand wildcards in given filepatterns.

    Arguments
    *args -- first arg: filepatterns as list or one single filepattern,
        second arg (optional): a function to combine wildcard values
        (itertools.product per default)
    **wildcards -- the wildcards as keyword arguments
        with their values as lists. If allow_missing=True is included
        wildcards in filepattern without values will stay unformatted.
    """

class Namedlist[T](list[T | Iterable[T]]): ...

class InputFiles[T](Namedlist[T]):
    def __getattribute__(self, __name: str) -> str: ...

class OutputFiles[T](Namedlist[T]):
    def __getattribute__(self, __name: str) -> str: ...

class Params[T](Namedlist[T]): ...

class Snakemake:
    input: InputFiles[str]
    output: OutputFiles[str]
    params: Params[str]

def main(argv: list[str] = ...) -> None: ...
def get_argument_parser(profiles: list[str] | None = ...) -> ArgumentParser: ...
