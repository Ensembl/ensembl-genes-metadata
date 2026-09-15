"""Sphinx configuration for ensembl-genes-metadata."""

import os
import sys

sys.path.insert(0, os.path.abspath("../../src/python"))

project = "ensembl-genes-metadata"
copyright = "2026, Genebuild"
author = "Genebuild"

extensions = [
    "autoapi.extension",
    "myst_parser",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
]

# sphinx-autoapi parses the source tree directly, so it doesn't need
# gb_metadata's runtime dependencies installed to build the API docs.
autoapi_type = "python"
autoapi_dirs = ["../../src/python"]
autoapi_root = "api"
autoapi_add_toctree_entry = True
autoapi_options = [
    "members",
    "undoc-members",
    "show-inheritance",
    "show-module-summary",
]

source_suffix = {
    ".rst": "restructuredtext",
    ".md": "markdown",
}

exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

html_theme = "sphinx_rtd_theme"
