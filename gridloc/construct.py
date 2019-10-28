from pathlib import Path
from numpy import ones, NaN, pi
from logging import getLogger

from .geometry import compute_neighbor
from .search import find_new_pos_1d, find_new_pos_2d
from .generators import index_up_down, index_spiral

lg = getLogger(__name__)

CWD = Path(__file__).parent
DATA_PATH = CWD.parent / 'tests' / 'data'


def construct_grid(surf, start_vert, n_rows, n_cols, rotation=0,
                   index='up_down'):

    radians = rotation / 180 * pi
    neighbors = compute_neighbor(n_rows, n_cols)

    grid = ones((n_rows, n_cols, 2, 3))
    grid.fill(NaN)

    if index == 'up_down':
        g = index_up_down(n_rows, n_cols)
    elif index == 'spiral':
        g = index_spiral(n_rows, n_cols)

    [x_start, y_start] = next(g)
    grid[x_start, y_start, 0, :] = surf['pos'][start_vert, :]
    grid[x_start, y_start, 1, :] = surf['pos_norm'][start_vert, :]

    for x, y in g:
        n_neighbors = (neighbors[x, y, :, 0] > -1).sum()
        if n_neighbors == 0:
            raise ValueError('It cannot have zero neighbors')
        elif n_neighbors == 1:
            grid[x, y, :, :] = find_new_pos_1d(x, y, grid, neighbors, surf, radians=radians)
        else:
            grid[x, y, :, :] = find_new_pos_2d(x, y, grid, neighbors, surf)

    return grid


if __name__ == '__main__':

    construct_grid(
        str(DATA_PATH / 'rh.smooth'),
        2288,
        16,
        8)
