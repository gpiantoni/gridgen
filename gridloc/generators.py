from numpy import where


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

    Returns
    -------
    generator of tuples
        each tuple contains the index for the row (x) and the column (y)

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
