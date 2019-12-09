from numpy import zeros, dtype
from numpy.linalg import norm
from scipy.stats import norm as normal_dist


def compute_vasculature(grid, angio):

    d_ = dtype([
        ('label', '<U256'),   # labels cannot be longer than 256 char
        ('vasculature', 'f4'),
        ])

    distance = zeros((grid.shape[0], grid.shape[1]), dtype=d_)
    for i_x in range(grid.shape[0]):
        for i_y in range(grid.shape[1]):
            weights = normal_dist.pdf(norm(angio['pos'] - grid['pos'][i_x, i_y], axis=1), scale=2)
            distance['vasculature'][i_x, i_y] = (weights * angio['value']).sum()

    distance['label'] = grid['label']
    return distance
