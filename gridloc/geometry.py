from numpy import ones

from .generators import index_up_down, index_spiral


def compute_neighbor(n_rows, n_cols, index='up_down'):
    """make sure that any indexing does not have more than 2 neighbors
    """

    neigh = -1 * ones((n_rows, n_cols, 2, 2), dtype=int)

    if index == 'up_down':
        g = index_up_down(n_rows, n_cols)
    elif index == 'spiral':
        g = index_spiral(n_rows, n_cols)

    [x_start, y_start] = next(g)
    neigh[x_start, y_start, 0, :] = [x_start, y_start]

    neighbors = [
        [1, 0],
        [-1, 0],
        [0, 1],
        [0, -1],
        ]
    for x, y in g:
        i = 0
        for n_x, n_y in neighbors:
            if (x + n_x) < 0 or (x + n_x) >= n_rows:
                continue
            if (y + n_y) < 0 or (y + n_y) >= n_cols:
                continue
            if neigh[x + n_x, y + n_y, 0, 0] == -1:
                continue
            neigh[x, y, i, :] = x + n_x, y + n_y
            i += 1

    neigh[x_start, y_start, 0, :] = -1

    return neigh
