"""Script to generate a Snakebids project."""

import argparse
import re
import sys
from pathlib import Path

import copier
import more_itertools as itx
from colorama import Fore, Style

import snakebids
from snakebids.app import SnakeBidsApp
from snakebids.cli import add_dynamic_args


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
    print(
        f"Creating Snakebids app at {Fore.GREEN}{output}{Fore.RESET}", file=sys.stderr
    )
    print(file=sys.stderr)
    try:
        copier.run_copy(
            str(Path(itx.first(snakebids.__path__), "project_template")),
            output,
            data={"app_full_name": output.name},
            unsafe=True,
        )
    except KeyboardInterrupt:
        print(f"{Fore.RED}Aborted!{Fore.RESET}", file=sys.stderr)
        sys.exit(1)


def create_descriptor(args: argparse.Namespace) -> None:
    """Implement the ``snakebids boutiques`` command."""
    app = SnakeBidsApp(args.app_dir.resolve())
    add_dynamic_args(app.parser, app.config["parse_args"], app.config["pybids_inputs"])
    app.create_descriptor(args.out_path)
    print(f"Boutiques descriptor created at {args.out_path}")


def gen_parser() -> argparse.ArgumentParser:
    """Generate the CLI parser for ``snakebids``."""
    parser = argparse.ArgumentParser(
        description="Perform administrative Snakebids tasks."
    )
    subparsers = parser.add_subparsers(required=True, dest="command")

    parser_create = subparsers.add_parser("create", help="Create a new Snakebids app.")
    parser_create.add_argument("output_dir", nargs="?", default=".")
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
