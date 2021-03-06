"""Function to build a grid from a starting point, in 2D
"""
from numpy import NaN, zeros, array, fliplr, flipud, vstack
from logging import getLogger
from .utils import DTYPE_MESH, WIRE


lg = getLogger(__name__)


def make_grid_with_labels(n_rows, n_columns, direction, chan_pattern='{}'):
    """Create a regular 2D grid, with labels and wires

    Parameters
    ----------
    n_rows : int
        number of rows
    n_columns : int
        number of columns
    direction : str
        specify direction of the labels. One of:
        'LRTB', 'RLTB', 'LRBT', 'RLBT', 'TBLR', 'TBRL', 'BTLR', 'BTRL'
        (where L -> Left, R -> Right, T -> Top, B -> Bottom)
    chan_pattern : str
        pattern to convert number into label string. Examples:

            'chan{}' -> 'chan1', 'chan2', 'chan3', ...
            'chan{:02d}' -> 'chan01', 'chan02', 'chan03', ...
            'elec{:03d}' -> 'elec001', 'elec002', 'elec003', ...

    Returns
    -------
    ndarray of shape (n_rows + 1, n_columns) with fields:
        - label : str (labels)
        - pos : 3 floats (specifying the x, y, z position)
        - norm : 3 floats (specifying the normals)
        - done : bool (where the position has been computed or not)

    The labels in the bottom row are all WIRE
    """
    grid = make_grid(n_rows, n_columns)

    if direction[0] in ('L', 'R'):
        order = 'C'
    else:
        order = 'F'

    labels = array([chan_pattern.format(x + 1) for x in range(n_rows * n_columns)]).reshape(n_rows, n_columns, order=order)

    if 'RL' in direction:
        labels = fliplr(labels)
    if 'BT' in direction:
        labels = flipud(labels)
    grid['label'] = labels

    wires = zeros((1, n_columns), dtype=DTYPE_MESH)
    wires['label'] = WIRE

    return vstack([grid, wires])


def make_grid(n_rows, n_columns):
    """Create an empty regular 2D grid

    Parameters
    ----------
    n_rows : int
        number of rows
    n_columns : int
        number of columns

    Returns
    -------
    ndarray of shape (n_rows, n_columns) with fields:
        - label : str (labels)
        - pos : 3 floats (specifying the x, y, z position)
        - norm : 3 floats (specifying the normals)
        - done : bool (where the position has been computed or not)
    """
    grid = zeros((n_rows, n_columns), dtype=DTYPE_MESH)
    grid['pos'].fill(NaN)
    grid['norm'].fill(NaN)
    grid['done'] = False

    return grid
