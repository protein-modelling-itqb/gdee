# Configuration file for the Sphinx documentation builder.

import os
import sys

# Add the package to the Python path
sys.path.insert(0, os.path.abspath('../../package'))

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = "GDEE"
copyright = "2025, Caio S. Souza"
author = "Caio S. Souza"
release = "1.0.0"

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx.ext.autosummary",
    "sphinx.ext.coverage",
    "sphinx.ext.mathjax",
    "sphinx.ext.githubpages",
    "sphinx.ext.viewcode",
    "sphinx.ext.intersphinx",
]

# Avoid RTD failures on optional / non-trivial native deps during autodoc imports
autodoc_mock_imports = [
    "modeller",
    "MDAnalysis",
    "mpi4py",
    "oddt",
    "openbabel",
]

templates_path = ["_templates"]
exclude_patterns = []


# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = "sphinx_rtd_theme"  # Better theme than alabaster
html_static_path = ["_static"]

# -- Extension configuration -------------------------------------------------
autodoc_default_options = {
    'members': True,
    'member-order': 'bysource',
    'special-members': '__init__',
    'undoc-members': True,
    'exclude-members': '__weakref__'
}

# Napoleon settings for Google/NumPy style docstrings
napoleon_google_docstring = True
napoleon_numpy_docstring = True
napoleon_include_init_with_doc = False
napoleon_include_private_with_doc = False

# Intersphinx mapping for cross-references
intersphinx_mapping = {
    'python': ('https://docs.python.org/3/', None),
    'numpy': ('https://numpy.org/doc/stable/', None),
    'mdanalysis': ('https://docs.mdanalysis.org/stable/', None),
    'biopython': ('https://biopython.org/docs/latest/api/', None),
}

# Autosummary settings
autosummary_generate = True
