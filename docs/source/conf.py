# -*- coding: utf-8 -*-
#
# Configuration file for the Sphinx documentation builder.
#
# This file does only contain a selection of the most common options. For a
# full list see the documentation:
# http://www.sphinx-doc.org/en/master/config

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
# import os
# import sys
# sys.path.insert(0, os.path.abspath('.'))
import time
import pathlib as pl

from pypath._metadata import __author__, __license__, __version__

# -- Project information -----------------------------------------------------

project = u'pypath'
copyright = ', '.join(
    [time.strftime('%Y'), 
    (", ").join(__author__)]
)
author = (", ").join(__author__)

# The short X.Y version
version = __version__
# The full version, including alpha/beta/rc tags
release = __version__

# thank you stupid sphinx, thank you stupid github :(((
readme_lines = []
readme = pl.Path().absolute().parents[1].joinpath('README.rst')

if readme.exists():

    with readme.open('r') as fp:

        readme_lines = fp.readlines()[4:]

with open('index.rst', 'w') as fp:

    fp.write('==================\nWelcome to pypath!\n==================\n\n')
    fp.write(''.join(readme_lines))

# -- General configuration ---------------------------------------------------

# If your documentation needs a minimal Sphinx version, state it here.
#
# needs_sphinx = '1.0'

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    'sphinx.ext.napoleon',
    'sphinx.ext.autosummary',
    'sphinx.ext.ifconfig',
    'sphinx.ext.viewcode',
    'sphinx.ext.githubpages',
    'sphinx.ext.mathjax',
    'sphinx_copybutton',
    'sphinx_last_updated_by_git',
    'sphinxcontrib.fulltoc',
    'sphinx_remove_toctrees',
    'nbsphinx',
    'IPython.sphinxext.ipython_console_highlighting',
]

autosummary_generate = True
remove_from_toctrees = ['api/*']

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# The suffix(es) of source filenames.
# You can specify multiple suffix as a list of string:
#
# source_suffix = ['.rst', '.md']
source_suffix = '.rst'

# The master toctree document.
master_doc = 'contents'

# The language for content autogenerated by Sphinx. Refer to documentation
# for a list of supported languages.
#
# This is also used if you do content translation via gettext catalogs.
# Usually you set "language" from the command line for these cases.
language = 'en'

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path .
exclude_patterns = []

# The name of the Pygments (syntax highlighting) style to use.
pygments_style = 'manni'


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_permalinks_icon = '§'
html_theme = 'insipid'
html_theme = 'pydata_sphinx_theme'
html_logo = '_static/img/omnipath-logo.svg'
html_favicon = '_static/img/omnipath-logo.svg'

# Theme options are theme-specific and customize the look and feel of a theme
# further.  For a list of options available for each theme, see the
# documentation.
#
# html_theme_options = {}

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']

html_css_files = [
    'css/pypathdocs.css',
]

html_theme_options = {
    'show_toc_level': 2,
    'icon_links': [
        {
            'name': 'PyPI',
            'url': 'https://pypi.org/project/pypath-omnipath',
            'icon': 'fa-solid fa-box',
        },
        {
            'name': 'OmniPath',
            'url': 'https://omnipathdb.org/',
            'icon': 'fa-solid fa-university',
        },
        {
            'name': 'Saez Lab',
            'url': 'https://saezlab.org/',
            'icon': 'fa-solid fa-users',
        },
    ],
    'logo': {
        'text': 'Pypath',
        'image_light': 'img/omnipath-logo.svg',
        'image_dark': 'img/omnipath-logo.svg',
        'alt_text': 'Pypath',
    },
    'show_toc_level': 2,
    'github_url': 'https://github.com/saezlab/pypath',
    'twitter_url': 'https://twitter.com/omnipathdb',
    'header_links_before_dropdown': 5,
    'navbar_align': 'left',
    'navbar_center': ['navbar-nav'],
}

# Custom sidebar templates, must be a dictionary that maps document names
# to template names.
#
# The default sidebars (for documents that don't match any pattern) are
# defined by theme itself.  Builtin themes are using these templates by
# default: ``['localtoc.html', 'relations.html', 'sourcelink.html',
# 'searchbox.html']``.
#
# html_sidebars = {}

html_sidebars = {
    'notebooks/*': [
        'localtoc.html',
    ]
}

# -- Options for HTMLHelp output ---------------------------------------------

# Output file base name for HTML help builder.
htmlhelp_basename = 'pypathdoc'


# -- Options for LaTeX output ------------------------------------------------

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
    (master_doc, 'pypath.tex', u'pypath Documentation',
     u'Dénes Türei', 'manual'),
]


# -- Options for manual page output ------------------------------------------

# One entry per manual page. List of tuples
# (source start file, name, description, authors, manual section).
man_pages = [
    (master_doc, 'pypath', u'pypath Documentation',
     [author], 1)
]


# -- Options for Texinfo output ----------------------------------------------

# Grouping the document tree into Texinfo files. List of tuples
# (source start file, target name, title, author,
#  dir menu entry, description, category)
texinfo_documents = [
    (master_doc, 'pypath', u'pypath Documentation',
     author, 'pypath', 'One line description of project.',
     'Miscellaneous'),
]


# -- Extension configuration -------------------------------------------------
