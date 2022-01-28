#!/usr/bin/env python3
from pathlib import Path

from snakebids.app import SnakeBidsApp


def get_parser():
    """Exposes parser for sphinx doc generation, cwd is the docs dir"""
    app = SnakeBidsApp.from_filesystem("../", skip_parse_args=True)
    return app.parser


def main():
    # to get repository root
    app = SnakeBidsApp.from_filesystem(Path(__file__).resolve().parents[1])
    app.run_snakemake()


if __name__ == "__main__":
    main()
