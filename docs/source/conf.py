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

# autodoc options
autosummary_generate = True

# todo settings
todo_include_todos = True
todo_link_only = True
