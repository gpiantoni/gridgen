from numpy import ones, array, zeros, argmin
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


def compute_neighbor(n_rows, n_cols, index='up_down'):
    """
    Deprecated

    make sure that any indexing does not have more than 2 neighbors
    Also, when you have two neighbors, the order is important. If you draw a
    line from the first neighbor to the second neighbor looking from above,
    the target position should be on the right
    """
    # include padding, to avoid indexerror
    neigh = -1 * ones((n_rows + 2, n_cols + 2, 2, 2), dtype=int)

    if index == 'up_down':
        g = index_up_down(n_rows, n_cols)
    elif index.startswith('corner'):
        g = index_corner(n_rows, n_cols, start=index[7:])
    else:
        raise NotImplementedError(f'index {index} not implemented')

    [x_start, y_start] = next(g)
    neigh[x_start + 1, y_start + 1, 0, :] = [x_start, y_start]

    neighbors = array([
        [-1, 0],
        [0, 1],
        [1, 0],
        [0, -1],
        ])

    for x, y in g:
        n = neighbors + (x, y)
        idx = (neigh[n[:, 0] + 1, n[:, 1] + 1, 0, 0] != -1)
        if idx.sum() == 1:
            neigh[x + 1, y + 1, 0, :] = n[idx]
        elif idx.sum() == 2:
            if (idx == array([True, False, False, True])).all():
                neigh[x + 1, y + 1, 0, :] = n[3]
                neigh[x + 1, y + 1, 1, :] = n[0]
            else:
                neigh[x + 1, y + 1, :, :] = n[idx]
        else:
            raise IndexError('There is something wrong')

    neigh[x_start + 1, y_start + 1, 0, :] = -1
    neigh = neigh[1:-1, 1:-1, :, :]  # remove padding

    return neigh
