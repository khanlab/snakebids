#!/usr/bin/env python3
from pathlib import Path

from snakebids.app import SnakeBidsApp
from snakebids.cli import create_parser


def get_parser():
    """Exposes parser for sphinx doc generation, cwd is the docs dir"""
    return create_parser()


def main():
    # to get repository root
    app = SnakeBidsApp.from_filesystem(Path(__file__).resolve().parents[1])
    app.run_snakemake()


if __name__ == "__main__":
    main()
