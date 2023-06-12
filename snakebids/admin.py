"""Script to generate a Snakebids project."""

import argparse
from pathlib import Path

from cookiecutter.main import cookiecutter  # type: ignore

import snakebids
from snakebids.app import SnakeBidsApp
from snakebids.cli import add_dynamic_args


def create_app(args: argparse.Namespace) -> None:
    cookiecutter(
        str(Path(list(snakebids.__path__)[0]) / "project_template"),
    )


def create_descriptor(args: argparse.Namespace) -> None:
    # pylint: disable=unsubscriptable-object
    app = SnakeBidsApp(args.app_dir.resolve())
    add_dynamic_args(app.parser, app.config["parse_args"], app.config["pybids_inputs"])
    app.create_descriptor(args.out_path)
    print(f"Boutiques descriptor created at {args.out_path}")


def gen_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Perform administrative Snakebids tasks."
    )
    subparsers = parser.add_subparsers(required=True, dest="command")

    parser_create = subparsers.add_parser("create", help="Create a new Snakebids app.")
    parser_create.add_argument("output_dir", nargs="?", default=".")
    parser_create.set_defaults(func=create_app)

    parser_boutiques = subparsers.add_parser(
        "boutiques", help="Create a Boutiques descriptor for an existing Snakebids app."
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
    """Invoke Cookiecutter on the Snakebids project template."""

    parser = gen_parser()
    args = parser.parse_args()
    args.func(args)
