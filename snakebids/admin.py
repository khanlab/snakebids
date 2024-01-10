"""Script to generate a Snakebids project."""

import argparse
import re
import sys
from pathlib import Path

import copier
import more_itertools as itx
from colorama import Fore, Style

import snakebids
from snakebids.utils.utils import text_fold


def create_app(args: argparse.Namespace) -> None:
    """Implement the ``snakebids create`` command."""
    output = Path(args.output_dir).resolve()
    if not output.parent.exists():
        print(
            f"{Fore.RED}{Style.BRIGHT}{output.parent}{Style.RESET_ALL}{Fore.RED} does "
            f"not exist{Fore.RESET}",
            file=sys.stderr,
        )
        sys.exit(1)
    if not re.match(r"^[a-zA-Z_][a-zA-Z_0-9]*$", output.name):
        print(
            f"{Fore.RED}Output directory name {Style.BRIGHT}{output.name}"
            f"{Style.RESET_ALL}{Fore.RED} is not a valid python module name",
            file=sys.stderr,
        )
        sys.exit(1)

    data = {"app_full_name": output.name}

    if args.snakebids_version is not None:
        version = args.snakebids_version
        if ("@" not in version and ";" in version) or (
            "@" in version and " ;" in version
        ):
            print(
                f"{Fore.RED}Snakebids version may not specify markers{Style.RESET_ALL}",
                file=sys.stderr,
            )
            sys.exit(1)
        if (
            "@" in version
            and version[1:].lstrip().startswith("git+")
            and "@" in version[1:]
        ):
            print(
                f"{Fore.RED}Credentials and rev specifiers in git requirement "
                f"specifications are not supported{Style.RESET_ALL}",
                file=sys.stderr,
            )
            sys.exit(1)
        if version.strip().startswith("["):
            print(
                f"{Fore.RED}Snakebids version may not specify extras{Style.RESET_ALL}",
                file=sys.stderr,
            )
            sys.exit(1)

        data["snakebids_version"] = args.snakebids_version

    print(
        f"Creating Snakebids app at {Fore.GREEN}{output}{Fore.RESET}", file=sys.stderr
    )
    print(file=sys.stderr)
    try:
        copier.run_copy(
            str(Path(itx.first(snakebids.__path__), "project_template")),
            output,
            data=data,
            unsafe=True,
        )
    except KeyboardInterrupt:
        print(f"{Fore.RED}Aborted!{Fore.RESET}", file=sys.stderr)
        sys.exit(1)


def create_descriptor(args: argparse.Namespace) -> None:
    """Implement the ``snakebids boutiques`` command."""
    print("Boutiques descriptor generation is temporarily disabled", file=sys.stderr)
    sys.exit(1)


def gen_parser() -> argparse.ArgumentParser:
    """Generate the CLI parser for ``snakebids``."""
    parser = argparse.ArgumentParser(
        description="Perform administrative Snakebids tasks."
    )
    subparsers = parser.add_subparsers(required=True, dest="command")

    parser_create = subparsers.add_parser("create", help="Create a new Snakebids app.")
    parser_create.add_argument("output_dir", nargs="?", default=".")
    parser_create.add_argument(
        "--snakebids-version",
        default=None,
        metavar="VERSION_SPECIFIER",
        help=text_fold(
            """
            Specify snakebids version requirement. Supports either a valid version
            specifier (e.g. `>=x.x.x`, `==a.b.c`) or a url prepended with `@` (e.g. `@
            https://...`). Paths can be specified with `@ file:///absolute/path/...`.
            Markers and extras may not be specified.
            """
        ),
    )
    parser_create.set_defaults(func=create_app)

    parser_boutiques = subparsers.add_parser(
        "boutiques",
        help="Create a Boutiques descriptor for an existing Snakebids app.",
    )
    parser_boutiques.add_argument(
        "out_path",
        help="Path for the output Boutiques descriptor. Should be a .json file.",
        type=Path,
    )
    parser_boutiques.add_argument(
        "--app_dir",
        help="Location of the Snakebids app. Defaults to the current directory.",
        type=Path,
        default=".",
    )
    parser_boutiques.set_defaults(func=create_descriptor)

    return parser


def main() -> None:
    """Invoke Snakebids cli."""
    parser = gen_parser()
    args = parser.parse_args()
    args.func(args)
