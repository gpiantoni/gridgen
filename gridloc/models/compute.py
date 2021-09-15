"""Functions to compute the actual fitting of the grid onto the brain surface
"""
from logging import getLogger

from .morphology import compute_morphology
from .functional import compute_functional

lg = getLogger(__name__)


def compute_model(grid, mri, morphology={}, functional={}):

    if mri['pial'] is None:
        morpho = None
    else:
        morpho = compute_morphology(grid, mri['pial'], **morphology)

    if mri['functional'] is None:
        func = None
    else:
        func = compute_functional(grid, mri['func'], **functional)

    model = {
        'grid': grid,
        'morphology': morpho,
        'functional': func,
        }

    return model
