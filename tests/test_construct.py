
from pathlib import Path
from gridloc.construct import construct_grid
from gridloc.io import read_surf, export_grid

from numpy.testing import assert_array_almost_equal
from numpy import array

TEST_PATH = Path(__file__).resolve().parent
DATA_PATH = TEST_PATH / 'data'


def test_construct():
    surf_file = DATA_PATH / 'lh_smooth.pial'
    grid_file = DATA_PATH / 'grid_020.fcsv'

    surf = read_surf(surf_file)
    grid = construct_grid(surf, 33154, 16, 8, rotation=20)
    assert_array_almost_equal(
        grid[0, 0, 0, :],
        array([-49.076, -11.536, 46.869]),
        decimal=3)

    export_grid(grid, grid_file, 'slicer')
    export_grid(grid, grid_file, 'freeview')
