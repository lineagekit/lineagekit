import os
import sys

sys.path.insert(0, os.path.abspath('../../src/core/'))

# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'lineagekit'
copyright = '2025, Andrii Serdiuk'
author = 'Andrii Serdiuk'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.napoleon',
    'sphinx_autodoc_typehints',
]

autodoc_default_options = {
    'members': True,
    'undoc-members': True,
    'private-members': False,  # Excludes private members by default
    'special-members': '__init__',
}

autodoc_typehints = "description"  # Show types in description instead of signature
typehints_defaults = "braces"      # Show default values nicely

napoleon_google_docstring = True
napoleon_numpy_docstring = False

templates_path = ['_templates']
exclude_patterns = []

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_book_theme'
html_static_path = ['_static']
