from logging import getLogger

from numpy import (power,
                   zeros,
                   NaN,
                   nansum,
                   dtype,
                   )
from numpy.linalg import norm
from scipy.stats import norm as normdistr
from ..io import WIRE

lg = getLogger(__name__)

# 1 sigma = 0.6065306597126334


def compute_functional(grid, func, distance=None, kernel=None):

    d_ = dtype([
        ('label', '<U256'),   # labels cannot be longer than 256 char
        ('value', 'f4'),
        ])

    dist = zeros((grid.shape[0], grid.shape[1]), dtype=d_)
    dist['label'] = grid['label']
    dist['value'].fill(NaN)

    for i_x in range(grid.shape[0]):
        for i_y in range(grid.shape[1]):
            if grid['label'][i_x, i_y] == WIRE:
                continue
            dist['value'][i_x, i_y] = compute_value_at_elec(
                grid['pos'][i_x, i_y], func, distance, kernel)

    return dist


def compute_value_at_elec(pos, func, distance='gaussian', kernel=8):
    """

    TODO
    ----
    How to normalize:
    You need to think hard if you need to normalize value (because we only
    include voxels > 0)
    """
    dist = norm(func['pos'] - pos, axis=1)

    if distance == 'gaussian':
        m = normdistr.pdf(dist, scale=kernel)

    elif distance == 'sphere':
        m = zeros(dist.shape)
        m[dist <= kernel] = 1

    elif distance == 'inverse':
        m = power(dist, -1 * kernel)

    v = func['value'] * m
    return nansum(v)
