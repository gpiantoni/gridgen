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


from numpy import arange, where

def _mirror(n_x, start_x):
    for i_x in range(n_x):
        i_x += start_x
        if i_x >= n_x:
            i_x =  n_x + (start_x - i_x)  - 1
        yield i_x

def test_():

    n_x = 6
    n_y = 3

    d = arange(n_x * n_y, dtype=int) + 1
    d = d.reshape(n_x, n_y, order='F')
    d

    start_elec = 1
    order = 'shorter'

    start_x, start_y = where(d == start_elec)
    start_x = start_x.item()
    start_y = start_y.item()

def index_order():
    swap = False
    if order == 'longer':
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
                    print(d[i_long, i_short])
                else:
                    print(d[i_short, i_long])

    elif order == 'shorter':
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
                    print(d[i_long, i_short])
                else:
                    print(d[i_short, i_long])
