# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'milopy'
copyright = '2022, Emma Dann'
author = 'Emma Dann'
release = '0.1.0'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    "myst_nb",
    # "sphinx.ext.viewcode",
    # "sphinx.ext.autosummary",
    'autoapi.extension'
]

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']
autoapi_dirs = ['../milopy']

# # Generate the API documentation when building
# autosummary_generate = True
# autodoc_member_order = "bysource"

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = "sphinx_rtd_theme"
html_static_path = ['_static']

# -- Exclude tests ----
# This is the expected signature of the handler for this event, cf doc


def autodoc_skip_member_handler(app, what, name, obj, skip, options):
    # Basic approach; you might want a regex instead
    return name.startswith("test_")

# Automatically called by sphinx at startup


def setup(app):
    # Connect the autodoc-skip-member event from apidoc to the callback
    app.connect('autodoc-skip-member', autodoc_skip_member_handler)
