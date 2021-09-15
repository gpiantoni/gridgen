"""Functions to compute the actual fitting of the grid onto the brain surface
"""
from scipy.optimize import brute, minimize
from scipy.stats import spearmanr
from multiprocessing import Pool
from numpy.linalg import norm
from numpy import arange, array, argmin, corrcoef, zeros, dtype, intersect1d, NaN, unravel_index
from logging import getLogger
from datetime import datetime
from json import dump

try:
    import mkl
except ImportError:
    mkl = None

"""
from .geometry import search_grid
from .morphology.distance import compute_distance
from .vascular.sphere import compute_vasculature
from .construct import construct_grid
from .io import read_mri, write_tsv, WIRE, export_grid
from .viz import plot_results, to_div, to_html, plot_electrodes
from .examine import measure_distances, measure_angles
from .utils import be_nice, match_labels, normalize
"""

lg = getLogger(__name__)


def corr_ecog_model(x0, ecog, grid3d, initial, mri, fit, final=False):
    """Main model to minimize

    Note that when final=False, cc should be minimized (the smaller the better)
    while when final=True, cc should be as large as possible (more intuitive
    reading)

    Parameters
    ----------
    ranges : None
        ignored, but easier to pass when working with arguments

    """
    grid = construct_grid(
        mri["dura"],
        start_vert,
        initial["label"],
        ecog['label'],
        grid3d,
        rotation=initial['rotation'] + rotation)

    if mri['pial'] is not None:
        morpho = compute_distance(grid, mri['pial'], fit['distance'], fit['maximum_distance'])
        morpho['value'] = morpho['value'] ** (-1 * fit['penalty'])

    if mri['angio'] is not None:
        vasc = compute_vasculature(grid, mri['angio'])
        e, m, v = match_labels(ecog, morpho, vasc)[1:]

    return model


def remove_wires(model):
    """We select rows with WIRE, so that we keep the 2d shape of the models
    """
    i_keep = model['grid']['label'] != WIRE
    i_keep = i_keep.all(axis=1)
    for field in ['ecog', 'grid', 'morphology', 'vasculature']:
        if model[field] is not None:
            model[field] = model[field][i_keep, :]

    return model


def find_vertex(dura, ras):
    init_ras = array(ras)
    vertex = argmin(norm(dura['pos'] - init_ras, axis=1))
    vert_dist = norm(init_ras - dura['pos'][vertex])

    lg.debug(f'Target RAS: {init_ras}, vertex #{vertex} RAS: {dura["pos"][vertex]} (distance = {vert_dist:0.3}mm)')

    return vertex
