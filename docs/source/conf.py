# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

import pathlib
import sys

sys.path.append(str(pathlib.Path(__file__).absolute().parent / "ext"))

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = "smores"
copyright = "2022, Austin Mroz, Lukas Turcani"
author = "Austin Mroz, Lukas Turcani"
release = "0.1.0"

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    "sphinx.ext.doctest",
    "sphinx.ext.napoleon",
    "sphinx.ext.autosummary",
    "sphinx.ext.intersphinx",
    "sphinx.ext.viewcode",
]

autosummary_imported_members = True

intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "numpy": ("https://numpy.org/doc/stable/", None),
}

autodoc_member_order = "groupwise"
autodoc_typehints = "description"
autodoc_mock_imports = ["psi4"]
autoclass_content = "both"

templates_path = ["_templates"]
exclude_patterns: list[str] = []


# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = "furo"
html_theme_options = {
    "light_css_variables": {
        "color-problematic": "#ee5151",
        "color-foreground-primary": "#ffffffcc",
        "color-foreground-secondary": "#9ca0a5",
        "color-foreground-muted": "#81868d",
        "color-foreground-border": "#666666",
        "color-background-primary": "#131416",
        "color-background-secondary": "#1a1c1e",
        "color-background-hover": "#1e2124ff",
        "color-background-hover--transparent": "#1e212400",
        "color-background-border": "#303335",
        "color-background-item": "#444",
        "color-announcement-background": "#000000dd",
        "color-announcement-text": "#eeebee",
        "color-brand-primary": "#2b8cee",
        "color-brand-content": "#368ce2",
        "color-highlighted-background": "#083563",
        "color-guilabel-background": "#08356380",
        "color-guilabel-border": "#13395f80",
        "color-api-keyword": "#9ca0a5",
        "color-highlight-on-target": "#333300",
        "color-admonition-background": "#18181a",
        "color-card-border": "#1a1c1e",
        "color-card-background": "#18181a",
        "color-card-marginals-background": "#1e2124ff",
        "color-inline-code-background": "#1a1c1e",
    }
}
html_static_path = ["_static"]


pygments_style = "custom_dracula.CustomDraculaStyle"
pygments_dark_style = "custom_dracula.CustomDraculaStyle"
