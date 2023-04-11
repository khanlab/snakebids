# pylint: disable=all
# flake8: noqa
# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

import datetime
import importlib.metadata as ilm

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import os
import sys

sys.path.insert(0, os.path.abspath("../"))

# -- Project information -----------------------------------------------------

project = "Snakebids"
release = ilm.version("snakebids")
copyright = f"{datetime.date.today().year}, Ali R. Khan"
author = "Ali R. Khan"


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    "sphinxarg.ext",
    "sphinx.ext.intersphinx",
    "sphinx.ext.napoleon",
    "sphinx.ext.autodoc",
    "sphinxcontrib.asciinema",
    "myst_parser",
    "sphinx_copybutton",
    "sphinx_reredirects",
]


myst_enable_extensions = [
    "attrs_block",
]
myst_enable_extensions = [
    "attrs_block",
    "attrs_inline",
    "tasklist",
    "deflist",
    "fieldlist",
]
myst_number_code_blocks = ["python", "yaml"]


napoleon_google_docstring = False
napoleon_numpy_docstring = True

autodoc_member_order = "bysource"
autodoc_typehints = "description"
autosummary_imported_members = True


# Add any paths that contain templates here, relative to this directory.
templates_path = ["_templates"]

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]


master_doc = "index"

intersphinx_mapping = {
    "pybids": ("https://bids-standard.github.io/pybids/", None),
    "python": ("https://docs.python.org/3", None),
    "snakemake": ("https://snakemake.readthedocs.io/en/stable/", None),
}


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = "furo"
html_theme_options = {
    "source_repository": "https://github.com/akhanf/snakebids",
    "source_branch": "main",
    "source_directory": "docs/",
}

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ["_static"]
# templates_path = ["_templates"]

redirects = {
    "migration/0.5_to_0.6.md": "migration/0.5_to_0.8.html",
}

sphinxcontrib_asciinema_defaults = {
    "preload": 1,
    "rows": 24,
    "speed": 3,
}
