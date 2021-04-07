"""Script to generate a Snakebids project."""

import os

from cookiecutter.main import cookiecutter

import snakebids


def main():
    """Invoke Cookiecutter on the Snakebids project template."""
    cookiecutter(os.path.join(snakebids.__path__[0], "project_template/"))
