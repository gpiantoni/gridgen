"""Functions to compute the grid and model on surface
"""
from logging import getLogger

from ..grid3d import construct_grid, find_vertex
from .morphology import compute_morphology
from .functional import compute_functional

lg = getLogger(__name__)


def make_grid3d_model(mris, grid2d, grid3d, initial, morphology={}, functional={}):
    """Create a 3d grid from an initial position (in RAS) and MRI scans and
    make the full morphology and functional models (if present).

    Parameters
    ----------
    mris : dict
        fields "dura" (necessary), "pial" (optional), "func" (optional)
    grid2d : (n_rows, n_cols) array
        array with the 2d labels
    grid3d : dict
        parameters to create 3d grid
    initial : dict
        parameters where to start creating the grid
    morphology : dict
        if present, parameters to compute morphology model
    functional : dict
        if present, parameters to compute functional model

    Returns
    -------
    dict
        - grid : 3d grid with positions and normals
        - morphology : morphology model for each electrode
        - functional : functional model for each electrode
    """
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
    """Mmake the full morphology and functional models (if present).

    Parameters
    ----------
    mris : dict
        fields "dura" (necessary), "pial" (optional), "func" (optional)
    grid3d : (n_rows, n_cols) array
        3d grid with positions and normals for each electrode
    morphology : dict
        if present, parameters to compute morphology model
    functional : dict
        if present, parameters to compute functional model

    Returns
    -------
    dict
        - grid : 3d grid with positions and normals (same as input)
        - morphology : morphology model for each electrode
        - functional : functional model for each electrode
    """

    if mri['pial'] is None or not morphology:
        morpho = None
    else:
        morpho = compute_morphology(grid, mri['pial'], **morphology)

    if mri['func'] is None or not functional:
        func = None
    else:
        func = compute_functional(grid, mri['func'], **functional)

    model = {
        'grid': grid,
        'morphology': morpho,
        'functional': func,
        }

    return model
