# type: ignore

try:
    from snakemake.cli import get_argument_parser, main
    from snakemake.common import configfile
    from snakemake.common.configfile import load_configfile
except ImportError:
    import snakemake.io as configfile
    from snakemake import get_argument_parser, main
    from snakemake.io import load_configfile

from snakemake.io import expand

# Handle different snakemake versions for regex function
try:
    from snakemake.io import regex_from_filepattern
except ImportError:
    from snakemake.io import regex as regex_from_filepattern

from snakemake.script import Snakemake

__all__ = [
    "Snakemake",
    "configfile",
    "expand",
    "get_argument_parser",
    "load_configfile",
    "main",
    "regex_from_filepattern",
]
