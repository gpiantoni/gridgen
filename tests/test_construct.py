
from pathlib import Path
from gridloc.construct import construct_grid
from gridloc.io import read_surf, export_grid_to_3dslicer

from numpy.testing import assert_array_almost_equal
from numpy import array

TEST_PATH = Path(__file__).resolve().parent
DATA_PATH = TEST_PATH / 'data'


def test_construct():
    surf_file = DATA_PATH / 'rh_smooth.pial'
    grid_file = DATA_PATH / 'grid_020.fscv'

    surf = read_surf(surf_file)

    grid = construct_grid(surf, 2288, 16, 8, rotation=20)
    assert_array_almost_equal(
        grid[0, 0, 0, :],
        array([12.742, -78.207, 49.987]),
        decimal=3)

    export_grid_to_3dslicer(grid, grid_file)
