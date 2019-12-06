from pathlib import Path
from numpy import NaN, pi, dtype, zeros, array
from logging import getLogger

from .geometry import count_neighbors
from .search import find_new_pos_0d, find_new_pos_1d, find_new_pos_2d
from .generators import index_order

lg = getLogger(__name__)

CWD = Path(__file__).parent
DATA_PATH = CWD.parent / 'tests' / 'data'


def construct_grid(surf, start_vert, start_label, labels, rotation=0):

    radians = rotation / 180 * pi
    n_rows, n_cols = labels.shape

    # make sure that grid is empty
    grid = make_grid(n_rows, n_cols, '{}')
    grid['label'] = labels
    grid['pos'].fill(0)
    grid['norm'].fill(0)
    grid['done'] = False

    if start_label not in grid['label']:
        raise ValueError(f'"{start_label}" is not one of the labels of the grid')
    g = index_order(grid, start_label)
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
                find_new_pos_1d(grid, neighbors, surf, x, y, (opposite_x, opposite_y))

            else:
                find_new_pos_0d(grid, neighbors, surf, x, y, radians=radians)

        elif n_neighbors == 2:
            find_new_pos_2d(grid, neighbors, surf, x, y)

        else:
            raise ValueError(f'Electrode {x:d}-{y:d} has {n_neighbors:d} but it can only have one or two neighbors')

    return grid


def make_grid(n_rows, n_columns, chan_pattern):

    d_ = dtype([
        ('label', '<U256'),   # labels cannot be longer than 256 char
        ('pos', 'f4', (3, )),
        ('norm', 'f4', (3, )),
        ('done', 'bool'),
        ])
    grid = zeros((n_rows, n_columns), dtype=d_)
    grid['label'] = array([chan_pattern.format(x + 1) for x in range(n_rows * n_columns)]).reshape(n_rows, n_columns, order='F')
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
