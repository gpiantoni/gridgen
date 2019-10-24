
from pathlib import Path
from gridloc.construct import construct_grid

from numpy.testing import assert_array_almost_equal
from numpy import array

TEST_PATH = Path(__file__).resolve().parent
DATA_PATH = TEST_PATH / 'data'

def test_constrct():
    surf_file = DATA_PATH / 'rh.smooth'
    grid = construct_grid(surf_file, 2288, 16, 8)
    assert_array_almost_equal(
        grid[0, 0, 0, :],
        array([5.481, -78.183, 46.655]),
        decimal=3)
