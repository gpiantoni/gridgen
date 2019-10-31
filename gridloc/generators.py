def index_up_down(n_x, n_y):
    """Create a generator to return x and y indices starting from the center.
    See documentation for examples.

    Parameters
    ----------
    n_x : int
        number of rows
    n_y : int
        number of columns

    Returns
    -------
    generator of tuples
        each tuple contains the index for the row (x) and the column (y)
    """
    y = n_y // 2

    y_sign = -1 if n_y % 2 else 1

    for y_i in range(n_y):

        y += y_sign * y_i
        y_sign = -y_sign

        x_sign = -1 if n_x % 2 else 1
        x = n_x // 2

        for x_i in range(n_x):
            x += x_sign * x_i
            x_sign = -x_sign
            yield (x, y)


def index_corner(n_x, n_y, start='NW'):
    """Create a generator to return x and y indices starting from one of the four corners.
    See documentation for examples.

    Parameters
    ----------
    n_x : int
        number of rows
    n_y : int
        number of columns
    start : str
        starting corner (one of NW, NE, SW, SE)

    Returns
    -------
    generator of tuples
        each tuple contains the index for the row (x) and the column (y)
    """
    start = start.upper()

    for i_x in range(n_x):
        for i_y in range(n_y):
            if start == 'NW':
                yield i_x, i_y
            elif start == 'NE':
                yield i_x, n_y - i_y - 1
            elif start == 'SW':
                yield n_x - i_x - 1, i_y
            elif start == 'SE':
                yield n_x - i_x - 1, n_y - i_y - 1
