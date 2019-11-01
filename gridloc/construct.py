from pathlib import Path
from numpy import ones, NaN, pi, dtype, zeros, array
from logging import getLogger

from .geometry import count_neighbors
from .search import find_new_pos_0d, find_new_pos_2d
from .generators import index_order

lg = getLogger(__name__)

CWD = Path(__file__).parent
DATA_PATH = CWD.parent / 'tests' / 'data'


def construct_grid(surf, start_vert, n_rows, n_cols, rotation=0):

    radians = rotation / 180 * pi

    grid = make_grid(n_rows, n_cols)

    g = index_order(grid, 'elec001')
    [x_start, y_start] = next(g)
    grid['pos'][x_start, y_start] = surf['pos'][start_vert, :]
    grid['norm'][x_start, y_start] = surf['pos_norm'][start_vert, :]
    grid['done'][x_start, y_start] = True

    for x, y in g:
        n_neighbors, neighbors = count_neighbors(grid, x, y)

        if n_neighbors == 0:
            raise ValueError(f'Electrode {x:d}-{y:d} cannot have zero neighbors')

        elif n_neighbors == 1:
            opposite_x, opposite_y = 2 * neighbors[0] - (x, y)

            if (0 <= opposite_x < n_rows) and (0 <= opposite_y < n_cols) and grid['done'][opposite_x, opposite_y]:
                grid = find_new_pos_1d(grid, x, y)  # to do

            else:
                grid = find_new_pos_0d(grid, neighbors, surf, x, y, radians=radians)

        elif n_neighbors == 2:
            grid = find_new_pos_2d(x, y, grid, neighbors, surf)

        else:
            raise ValueError(f'Electrode {x:d}-{y:d} has {n_neighbors:d} but it can only have one or two neighbors')

    return grid


def make_grid(n_rows, n_cols):

    d_ = dtype([
        ('label', '<U256'),   # labels cannot be longer than 256 char
        ('pos', 'f4', (3, )),
        ('norm', 'f4', (3, )),
        ('done', 'bool'),
        ])
    grid = zeros((n_rows, n_cols), dtype=d_)
    grid['label'] = array([f'elec{x + 1:03d}' for x in range(n_rows * n_cols)]).reshape(n_rows, n_cols, order='F')
    grid['pos'].fill(NaN)
    grid['norm'].fill(NaN)
    grid['done'] = False
    return grid


if __name__ == '__main__':

    construct_grid(
        str(DATA_PATH / 'rh.smooth'),
        2288,
        16,
        8)
