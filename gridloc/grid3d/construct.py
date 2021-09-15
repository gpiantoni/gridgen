"""Function to build a grid from a starting point, in 2D or 3D
"""
from numpy import NaN, pi, dtype, zeros, array, fliplr, flipud, where, arange
from logging import getLogger

from .geometry import count_neighbors
from .search import find_new_pos_0d, find_new_pos_1d, find_new_pos_2d

lg = getLogger(__name__)


def construct_grid(surf, start_vert, start_label, labels, grid3d, rotation=0):
    """Construct 3D grid, based on a starting vertex on a surface

    Parameters
    ----------
    surf : dict
        surface with normals. It's better to use the smooth surface to have
        reasonable results
    start_vert : int
        index of the initial vertex
    start_label : str
        label of the electrode, which is assumed to be in the initial vertex
    labels : 2d array
        a 2d array with the channel labels
    grid3d : dict
        - interelec_distance : float
        - maximum_angle : float
    rotation : float
        rotation from the inferior-superior axis (clockwise, in degrees)

    Returns
    -------
    ndarray of shape (n_rows, n_columns) with fields:
        - label : str (labels)
        - pos : 3 floats (specifying the x, y, z position)
        - norm : 3 floats (specifying the normals)
        - done : bool (where the position has been computed or not)
    """
    grid3d['angles'] = arange(
        -1 * grid3d['maximum_angle'],
        grid3d['maximum_angle'],
        grid3d['step_angle'])

    radians = rotation / 180 * pi
    n_rows, n_cols = labels.shape

    # make sure that grid is empty
    grid = make_grid(n_rows, n_cols)
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
                find_new_pos_1d(grid, neighbors, surf, grid3d, x, y, (opposite_x, opposite_y))

            else:
                find_new_pos_0d(grid, neighbors, surf, grid3d, x, y, radians=radians)

        elif n_neighbors == 2:
            find_new_pos_2d(grid, neighbors, surf, grid3d, x, y)

        else:
            raise ValueError(f'Electrode {x:d}-{y:d} has {n_neighbors:d} but it can only have one or two neighbors')

    return grid


def index_order(indices, start_point, order='major'):
    """Create a generator to return x and y indices starting from the center.
    See documentation for examples.

    Parameters
    ----------
    d : 2d array
        array representing the grid, where each number in the array represents
        the electrode number
    start_point : int
        number of the electrode number used as starting point
    order : 'major' or 'minor'
        whether to start from the major (longer) axis or from the minor (shorter)
        axis

    Yields
    ------
    int
        index for the row (x)
    int
        index for the row (x)

    Notes
    -----
    In theory, we could start by going up or down (now it only goes down), but
    I don't see the need at this moment.
    """
    n_x, n_y = indices.shape
    start_x, start_y = where(indices['label'] == start_point)
    start_x = start_x.item()
    start_y = start_y.item()

    if order == 'major':
        if (n_x >= n_y):
            swap = False
            n_long, n_short = n_x, n_y
            start_long, start_short = start_x, start_y
        else:
            swap = True
            n_long, n_short = n_y, n_x
            start_long, start_short = start_y, start_x

        for i_short in _mirror(n_short, start_short):
            for i_long in _mirror(n_long, start_long):
                if not swap:
                    yield i_long, i_short
                else:
                    yield i_short, i_long

    elif order == 'minor':
        if (n_x >= n_y):
            swap = True
            n_long, n_short = n_x, n_y
            start_long, start_short = start_x, start_y
        else:
            swap = False
            n_long, n_short = n_y, n_x
            start_long, start_short = start_y, start_x

        for i_long in _mirror(n_long, start_long):
            for i_short in _mirror(n_short, start_short):
                if swap:
                    yield i_long, i_short
                else:
                    yield i_short, i_long


def _mirror(n_x, start_x):
    for i_x in range(n_x):
        i_x += start_x
        if i_x >= n_x:
            i_x = n_x + (start_x - i_x) - 1
        yield i_x
