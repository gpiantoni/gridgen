from functools import partial
from itertools import product
from logging import getLogger
from multiprocessing import Pool
import warnings

from numpy import (ndindex,
                   array,
                   sum,
                   power,
                   zeros,
                   NaN,
                   isfinite,
                   nansum,
                   isnan,
                   dtype,
                   where,
                   )
from numpy.linalg import norm, inv
from scipy.stats import norm as normdistr
from nibabel.affines import apply_affine
from nibabel import load as nload
from ..io import WIRE

lg = getLogger(__name__)

# 1 sigma = 0.6065306597126334


def compute_functional(grid, func, distance=None, kernel=None):

    d_ = dtype([
        ('label', '<U256'),   # labels cannot be longer than 256 char
        ('value', 'f4'),
        ])

    distance = zeros((grid.shape[0], grid.shape[1]), dtype=d_)
    distance['label'] = grid['label']

    for i_x in range(grid.shape[0]):
        for i_y in range(grid.shape[1]):
            if grid['label'][i_x, i_y] == WIRE:
                continue

    distance = zeros((grid.shape[0], grid.shape[1]), dtype=d_)
    distance['value'].fill(NaN)

    for i_x in range(grid.shape[0]):
        for i_y in range(grid.shape[1]):
            if grid['label'][i_x, i_y] == WIRE:
                continue
            distance['value'][i_x, i_y] = compute_value_at_elec(
                grid['pos'][i_x, i_y], func, distance, kernel)

    return distance


def compute_value_at_elec(pos, func, distance='gaussian', kernel=8):
    """

    TODO
    ----
    How to normalize:
    You need to think hard if you need to normalize value (because we only
    include voxels > 0
    """
    distance = norm(func['pos'] - pos, axis=1)

    if True or distance == 'gaussian':
        m = normdistr.pdf(distance, scale=kernel)

    elif distance == 'sphere':
        m = zeros(dist_chan.shape)
        m[dist_chan <= KERNEL] = 1

    elif distance == 'inverse':
        m = power(dist_chan, -1 * KERNEL)

    v = func['value'] * m
    return nansum(v)
