from numpy import zeros, savetxt
from gridloc.generators import index_corner, index_up_down

from .paths import GENERATORS_PATH


def notest_generators_corner():

    n_rows = 2
    n_cols = 4

    for direction in ('nw', 'ne', 'sw', 'se'):

        d = zeros((n_rows, n_cols), dtype='int')
        for i, xy in enumerate(index_corner(n_rows, n_cols, direction)):
            d[xy] = i

        savetxt(
            GENERATORS_PATH / f'corner_{direction}.csv',
            d,
            delimiter=',',
            fmt='%d')


def notest_generators_updown():

    n_rows = 5
    n_cols = 4

    d = zeros((n_rows, n_cols), dtype='int')
    for i, xy in enumerate(index_up_down(n_rows, n_cols)):
        d[xy] = i

    savetxt(
        GENERATORS_PATH / f'center_updown.csv',
        d,
        delimiter=',',
        fmt='%d')
