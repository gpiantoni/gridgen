#!/usr/bin/env python3

project = 'gridloc'

html_theme = "pyramid"

extensions = [
    'sphinx.ext.autosummary',
    'sphinx.ext.autodoc',
    'sphinx.ext.githubpages',
    'sphinx.ext.napoleon',
    'sphinx.ext.todo',
    'sphinx.ext.viewcode',
    ]

# cleaner index for modindex
modindex_common_prefix = ['gridloc.', ]

# autodoc options
autosummary_generate = True

# napoleon options
napoleon_use_rtype = False
napoleon_use_param = False
napoleon_use_keyword = True

# todo options
todo_include_todos = True
todo_link_only = True
