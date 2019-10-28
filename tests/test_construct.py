
from pathlib import Path
from gridloc.construct import construct_grid
from gridloc.io import read_surf

from numpy.testing import assert_array_almost_equal
from numpy import array

TEST_PATH = Path(__file__).resolve().parent
DATA_PATH = TEST_PATH / 'data'


def test_construct():
    surf_file = DATA_PATH / 'rh_smooth.pial'
    surf = read_surf(surf_file)

    grid = construct_grid(surf, 2288, 16, 8)
    assert_array_almost_equal(
        grid[0, 0, 0, :],
        array([5.488, -78.169, 46.664]),
        decimal=3)
