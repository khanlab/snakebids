import argparse
from pathlib import Path

from snakebids.app import SnakeBidsApp
from snakebids.cli import add_dynamic_args


def gen_boutiques_parser():
    parser_boutiques = argparse.ArgumentParser(
        description="Generate a boutiques descriptor for this Snakebids app."
    )
    parser_boutiques.add_argument(
        "out_path", help="Path for the boutiques descriptor JSON file."
    )
    return parser_boutiques


def main():
    parser_boutiques = gen_boutiques_parser()
    args = parser_boutiques.parse_args()

    app = SnakeBidsApp.from_filesystem(Path(__file__).resolve().parents[1])
    add_dynamic_args(app.parser, app.config["parse_args"], app.config["pybids_inputs"])
    app.create_descriptor(args.out_path)
    print(f"Boutiques descriptor created at {args.out_path}")


if __name__ == "__main__":
    main()
