"""Functions to visualize 2d activity (ecog, models) and 3d activity (on surfaces)
"""
from .viz2d import plot_grid2d
from .results import plot_fitting, plot_grid3d
from .utils import to_html, to_div

__all__ = [
    'plot_grid2d',
    'plot_grid3d',
    'plot_fitting',
    'to_div',
    'to_html',
    ]
