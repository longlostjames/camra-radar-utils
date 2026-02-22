# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

import os
import sys

# Add the parent directory so autodoc can import the modules
sys.path.insert(0, os.path.abspath('..'))

# -- Project information -----------------------------------------------------
project   = 'CAMRa Radar Utils'
copyright = '2026, Chris Walden, UK Research & Innovation / National Centre for Atmospheric Science'
author    = 'Chris Walden'
release   = '1.3.0'

# -- General configuration ---------------------------------------------------
extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.autosummary',
    'sphinx.ext.napoleon',      # Google / NumPy docstring styles
    'sphinx.ext.viewcode',      # [source] links
    'sphinx.ext.intersphinx',   # cross-links to Python, NumPy etc.
    'sphinx.ext.todo',
    'sphinx_rtd_theme',
]

autosummary_generate = True

napoleon_google_docstring = True
napoleon_numpy_docstring  = True
napoleon_include_init_with_doc = False

intersphinx_mapping = {
    'python': ('https://docs.python.org/3', None),
    'numpy' : ('https://numpy.org/doc/stable/', None),
}

templates_path   = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

# Mock imports so autodoc can parse modules without needing all dependencies
autodoc_mock_imports = [
    'pyart',
    'netCDF4',
    'camra_utils',
]

# -- HTML output -------------------------------------------------------------
html_theme = 'sphinx_rtd_theme'

html_theme_options = {
    'navigation_depth': 4,
    'titles_only'     : False,
    'logo_only'       : False,
}

html_static_path = ['_static']

# -- Autodoc defaults --------------------------------------------------------
autodoc_default_options = {
    'members'          : True,
    'member-order'     : 'bysource',
    'special-members'  : '__init__',
    'undoc-members'    : False,
    'show-inheritance' : True,
}

autodoc_typehints = 'description'

todo_include_todos = True
