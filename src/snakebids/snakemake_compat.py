# type: ignore

try:
    from snakemake.cli import get_argument_parser, main
    from snakemake.common import configfile
    from snakemake.common.configfile import load_configfile
except ImportError:
    import snakemake.io as configfile
    from snakemake import get_argument_parser, main
    from snakemake.io import load_configfile

from snakemake.exceptions import WildcardError
from snakemake.io import expand
from snakemake.script import Snakemake

__all__ = [
    "load_configfile",
    "main",
    "expand",
    "Snakemake",
    "WildcardError",
    "configfile",
    "get_argument_parser",
]
