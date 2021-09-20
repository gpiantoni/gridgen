"""Functions to compute the grid and model on surface
"""
from logging import getLogger

from ..grid3d import construct_grid, find_vertex
from .morphology import compute_morphology
from .functional import compute_functional
from ..viz import plot_electrodes, to_div, to_html

lg = getLogger(__name__)


def make_grid3d_model(output, grid2d, mris, grid3d, initial, morphology={}, functional={}):

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

    fig = plot_electrodes(mris['pial'], grid, ref_label=initial['label'], angio=mris['angio'])
    grid_file = output / 'start_pos.html'
    to_html([to_div(fig), ], grid_file)


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
