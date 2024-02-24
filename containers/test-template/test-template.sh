#!/bin/sh

set -eu

method="$1"
script_name="$2"

cp -r /app/* /work
script="'${script_name}' tests/data tests/result participant -c1 --skip-bids-validation"
case "$method" in
    "setuptools" )
        python -m venv .venv
        .venv/bin/python -m pip install --no-color .
        PATH=".venv/bin:$PATH" eval "$script"
        ;;
    "poetry" )
        poetry install --no-ansi
        eval "poetry run $script"
        ;;
    "hatch" )
        hatch env create
        eval "hatch env run -- $script"
        ;;
    "pdm" )
        pdm install
        eval "pdm run $script"
        ;;
    "docs" )
        python -m venv .venv
        .venv/bin/python -m pip install .
        .venv/bin/python -m pip install -r docs/requirements.txt
        .venv/bin/sphinx-build docs build/docs -W
        ;;
    * )
        >&2 echo "Invalid method"
        exit 1
        ;;
esac
