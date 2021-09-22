"""Functions which are a copy of the matlab implementation or are necessary
for those matlab functions.
"""
from .compat import projectElectrodes
from .comparison import compare_to_matlab

__all__ = [
    'projectElectrodes',
    'compare_to_matlab',
    ]
