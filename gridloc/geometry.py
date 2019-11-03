from numpy import array, zeros, argmin
from numpy.linalg import norm
from gridloc.algebra import calc_plane_to_axis

NEIGHBORS = array([
    [-1, 0],
    [0, 1],
    [1, 0],
    [0, -1],
    ])


def count_neighbors(grid, x, y):
    n_rows, n_cols = grid.shape

    idx = zeros(4, bool)

    i = 0
    for n_x, n_y in NEIGHBORS:
        if (0 <= (x + n_x) < n_rows) and (0 <= (y + n_y) < n_cols) and grid['done'][x + n_x, y + n_y]:
            idx[i] = True
        i += 1

    n = idx.sum()

    """make sure that any indexing does not have more than 2 neighbors
    Also, when you have two neighbors, the order is important. If you draw a
    line from the first neighbor to the second neighbor looking from above,
    the target position should be on the right"""
    if (idx == array([True, False, False, True])).all():
        coords = NEIGHBORS[(3, 0), :] + (x, y)
    else:
        coords = NEIGHBORS[idx] + (x, y)
    return n, coords


def search_grid(surf, ref_vert, x, y):
    pos = surf['pos'][ref_vert, :]
    normal = surf['pos_norm'][ref_vert, :]

    coords_2d = array([x, y])
    plane = calc_plane_to_axis(normal)
    target = coords_2d @ plane + pos
    new_vector = argmin(norm(surf['pos'] - target, axis=1))

    return new_vector
