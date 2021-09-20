"""Functions to compute the grid and model on surface
"""
from logging import getLogger

from ..grid3d import construct_grid, find_vertex
from .morphology import compute_morphology
from .functional import compute_functional

lg = getLogger(__name__)


def make_grid3d_model(mris, grid2d, grid3d, initial, morphology={}, functional={}):

    vertex = find_vertex(mris['dura'], initial['RAS'])
    grid = construct_grid(
        mris["dura"],
        vertex,
        initial["label"],
        grid2d['label'],
        grid3d,
        rotation=initial['rotation'])

    model = compute_model(mris, grid, morphology, functional)

    return model


def compute_model(mri, grid, morphology={}, functional={}):

    if mri['pial'] is None or morphology is None:
        morpho = None
    else:
        morpho = compute_morphology(grid, mri['pial'], **morphology)

    if mri['func'] is None or functional is None:
        func = None
    else:
        func = compute_functional(grid, mri['func'], **functional)

    model = {
        'grid': grid,
        'morphology': morpho,
        'functional': func,
        }

    return model
