#!/bin/sh

set -eu

method="$1"
script_name="$2"

cp -r /app/* /work
script="'${script_name}' tests/data tests/result participant -c1 --skip-bids-validation"
case "$method" in
    "setuptools" )
        python -m venv .venv
        .venv/bin/python -m pip install .
        PATH=".venv/bin:$PATH" eval "$script"
        ;;
    "poetry" )
        poetry install
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
    * )
        >&2 echo "Invalid method"
        exit 1
        ;;
esac
