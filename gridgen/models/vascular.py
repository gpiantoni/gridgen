from numpy import zeros, dtype, NaN
from numpy.linalg import norm
from scipy.stats import norm as normal_dist

"""
from ..io import WIRE
"""


def compute_vasculature(grid, angio):
    """The lower the value (negative), the more suppression due to blood vessels
    """
    d_ = dtype([
        ('label', '<U256'),   # labels cannot be longer than 256 char
        ('value', 'f4'),
        ])

    distance = zeros((grid.shape[0], grid.shape[1]), dtype=d_)
    distance['value'].fill(NaN)
    for i_x in range(grid.shape[0]):
        for i_y in range(grid.shape[1]):
            if grid['label'][i_x, i_y] != WIRE:
                weights = normal_dist.pdf(norm(angio['pos'] - grid['pos'][i_x, i_y], axis=1), scale=2)
                distance['value'][i_x, i_y] = -1 * (weights * angio['value']).sum()

    distance['label'] = grid['label']
    return distance
