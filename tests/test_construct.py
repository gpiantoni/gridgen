from gridloc.construct import construct_grid, make_grid
from gridloc.io import read_surf, export_grid

from numpy.testing import assert_array_almost_equal
from numpy import array

from .paths import DATA_PATH, SMOOTH_FILE


def test_construct():
    grid_file = DATA_PATH / 'grid_020.fcsv'

    surf = read_surf(SMOOTH_FILE)
    grid2d = make_grid(16, 8)
    grid = construct_grid(surf, 33154, 'elec016', grid2d, rotation=20)
    assert_array_almost_equal(
        grid['pos'][0, 0],
        array([-28.26, -1.373, 54.757]),
        decimal=3)

    export_grid(grid, grid_file, 'slicer')
    export_grid(grid, grid_file, 'freeview')
