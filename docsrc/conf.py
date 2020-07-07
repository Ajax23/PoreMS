#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# PoreMS documentation build configuration file, created by
# sphinx-quickstart on Wed Apr 11 09:44:56 2018.
#
# This file is execfile()d with the current directory set to its
# containing dir.
#
# Note that not all possible configuration values are present in this
# autogenerated file.
#
# All configuration values have a default; values that are commented out
# serve to show the default.

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
# import sys
# import numpy as np
# import matplotlib as mpl
# import matplotlib.pyplot as plt
import sphinx_bootstrap_theme

# sys.path.append(os.path.abspath('../'))

# # Install package
# if sys.platform == "win32":
#     # os.system("pip install ../.")
#     os.system("..\tests\venv\Scripts\activate")
# else:
#     # os.system("pip install ../. &> /dev/null")
#     os.system("source ../tests/venv/bin/activate")


# -- General configuration ------------------------------------------------

# If your documentation needs a minimal Sphinx version, state it here.
#
# needs_sphinx = '1.0'

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = ['sphinx.ext.autodoc',
              'sphinx.ext.doctest',
              'sphinx.ext.todo',
              'sphinx.ext.mathjax',
              'sphinx.ext.autosummary',
              'numpydoc']

# Generate the API documentation when building
autosummary_generate = True
numpydoc_show_class_members = False

# autodoc_default_flags = ['members','private-members']

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# The suffix(es) of source filenames.
# You can specify multiple suffix as a list of string:
#
# source_suffix = ['.rst', '.md']
source_suffix = '.rst'

# The master toctree document.
master_doc = 'index'

# General information about the project.
project = 'PoreMS'
copyright = '2020, Hamzeh Kraus'
author = 'Hamzeh Kraus'

# The version info for the project you're documenting, acts as replacement for
# |version| and |release|, also used in various other places throughout the
# built documents.
#
# The short X.Y version.
version = '0.2'
# The full version, including alpha/beta/rc tags.
release = '0.2.0'

# The language for content autogenerated by Sphinx. Refer to documentation
# for a list of supported languages.
#
# This is also used if you do content translation via gettext catalogs.
# Usually you set "language" from the command line for these cases.
language = None

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This patterns also effect to html_static_path and html_extra_path
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store', '**.ipynb_checkpoints']

# The name of the Pygments (syntax highlighting) style to use.
pygments_style = 'sphinx'

# If true, `todo` and `todoList` produce output, else they produce nothing.
todo_include_todos = True


# -- Options for HTML output ----------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#

# html_logo = "pics/logo.svg"

html_theme = 'bootstrap'

html_theme_options = {
    'navbar_sidebarrel': False,
    'navbar_pagenav': True,
    'navbar_pagenav_name': "Page",
    'globaltoc_includehidden': "true",
    'navbar_class': "navbar",
    'navbar_fixed_top': "true",
    'source_link_position': "",
    'bootswatch_theme': "paper",
    'bootstrap_version': "3",
    'navbar_links': [("API", "api"), ("Process", "process"), ("Molecule", "molecule"), ("Pore", "pore"), ("Workflow", "workflow"),],
}

html_theme_path = sphinx_bootstrap_theme.get_html_theme_path()


# Theme options are theme-specific and customize the look and feel of a theme
# further.  For a list of options available for each theme, see the
# documentation.
#
# html_theme_options = {}

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']

# Custom sidebar templates, must be a dictionary that maps document names
# to template names.
#
# This is required for the alabaster theme
# refs: http://alabaster.readthedocs.io/en/latest/installation.html#sidebars
# html_sidebars = {
#     '**': [
#         'about.html',
#         'navigation.html',
#         'relations.html',  # needs 'show_related': True theme option to display
#         'searchbox.html',
#         'donate.html',
#     ]
# }


# -- Options for HTMLHelp output ------------------------------------------

# Output file base name for HTML help builder.
htmlhelp_basename = 'PoreMSdoc'


# -- Options for LaTeX output ---------------------------------------------

latex_elements = {
    # The paper size ('letterpaper' or 'a4paper').
    #
    # 'papersize': 'letterpaper',

    # The font size ('10pt', '11pt' or '12pt').
    #
    # 'pointsize': '10pt',

    # Additional stuff for the LaTeX preamble.
    #
    # 'preamble': '',

    # Latex figure (float) alignment
    #
    # 'figure_align': 'htbp',
}

# Grouping the document tree into LaTeX files. List of tuples
# (source start file, target name, title,
#  author, documentclass [howto, manual, or own class]).
latex_documents = [
    (master_doc, 'PoreMS.tex', 'PoreMS Documentation', 'Hamzeh Kraus', 'manual'),
]


# -- Options for manual page output ---------------------------------------

# One entry per manual page. List of tuples
# (source start file, name, description, authors, manual section).
man_pages = [
    (master_doc, 'porems', 'PoreMS Documentation', [author], 1)
]


# -- Options for Texinfo output -------------------------------------------

# Grouping the document tree into Texinfo files. List of tuples
# (source start file, target name, title, author,
#  dir menu entry, description, category)
texinfo_documents = [
    (master_doc, 'PoreMS', 'PoreMS Documentation',
     author, 'MS', 'One line description of project.',
     'Miscellaneous'),
]

# Add the 'copybutton' javascript, to hide/show the prompt in code
# examples, originally taken from scikit-learn's doc/conf.py


def setup(app):
    app.add_stylesheet('style.css')
