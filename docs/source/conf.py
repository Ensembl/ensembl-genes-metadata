"""Sphinx configuration for ensembl-genes-metadata."""

import os
import sys
from datetime import date

sys.path.insert(0, os.path.abspath("../../src/python"))

project = "ensembl-genes-metadata"
author = "Genebuild"
copyright = f"{date.today().year}, Genebuild"

extensions = [
    "autoapi.extension",
    "myst_parser",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
    "sphinx.ext.extlinks",
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

myst_enable_extensions = [
    "colon_fence",
    "deflist",
    "fieldlist",
]
myst_heading_anchors = 3

html_theme = "sphinx_rtd_theme"
html_title = "ensembl-genes-metadata"
html_logo = "../../docs_mkdocs/img/ebang.png"
html_theme_options = {
    "collapse_navigation": False,
    "navigation_depth": 4,
    "logo_only": True,
}

extlinks = {
    "repo": ("https://github.com/Ensembl/ensembl-genes-metadata/blob/main/%s", "%s"),
}

linkcheck_ignore = [
    r"https://github.com/Ensembl/ensembl-genes-metadata/.*",
]
