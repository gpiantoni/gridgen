from numpy import ones, array

from .generators import index_up_down


def compute_neighbor(n_rows, n_cols, index='up_down'):
    """make sure that any indexing does not have more than 2 neighbors
    Also, when you have two neighbors, the order is important. If you draw a
    line from the first neighbor to the second neighbor looking from above,
    the target position should be on the right
    """
    # include padding, to avoid indexerror
    neigh = -1 * ones((n_rows + 2, n_cols + 2, 2, 2), dtype=int)

    if index == 'up_down':
        g = index_up_down(n_rows, n_cols)
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
