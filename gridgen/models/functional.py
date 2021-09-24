from logging import getLogger

from numpy import (power,
                   zeros,
                   NaN,
                   nansum,
                   )
from numpy.linalg import norm
from scipy.stats import norm as normdistr

from ..utils import DTYPE, WIRE

lg = getLogger(__name__)


def compute_functional(grid, func, metric=None, kernel=None):

    dist = zeros((grid.shape[0], grid.shape[1]), dtype=DTYPE)
    dist['label'] = grid['label']
    dist['value'].fill(NaN)

    for i_x in range(grid.shape[0]):
        for i_y in range(grid.shape[1]):
            if grid['label'][i_x, i_y] == WIRE:
                continue
            dist['value'][i_x, i_y] = compute_value_at_elec(
                grid['pos'][i_x, i_y], func, metric, kernel)

    return dist


def compute_value_at_elec(pos, func, metric='gaussian', kernel=8):
    """

    TODO
    ----
    How to normalize:
    You need to think hard if you need to normalize value (because we only
    include voxels > 0)
    """
    dist = norm(func['pos'] - pos, axis=1)

    if metric == 'gaussian':
        m = normdistr.pdf(dist, scale=kernel)

    elif metric == 'sphere':
        m = zeros(dist.shape)
        m[dist <= kernel] = 1

    elif metric == 'inverse':
        m = power(dist, -1 * kernel)

    v = func['value'] * m
    return nansum(v)
