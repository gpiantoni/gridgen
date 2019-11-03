
from gridloc.construct import construct_grid
from gridloc.io import read_surf, export_grid

from numpy.testing import assert_array_almost_equal
from numpy import array

from .paths import DATA_PATH, SMOOTH_FILE


def test_construct():
    grid_file = DATA_PATH / 'grid_020.fcsv'

    surf = read_surf(SMOOTH_FILE)
    grid = construct_grid(surf, 33154, 16, 8, rotation=20)
    assert_array_almost_equal(
        grid[0, 0, 0, :],
        array([-10.036, -11.673, 63.269]),
        decimal=3)

    export_grid(grid, grid_file, 'slicer')
    export_grid(grid, grid_file, 'freeview')
