from numpy import array
from numpy.testing import assert_array_almost_equal

from gridloc.construct import construct_grid
from gridloc.fitting import fitting_brute, fitting_hop
from gridloc.io import read_surf
from gridloc.morphology.distance import compute_distance

from .paths import SMOOTH_FILE, PIAL_FILE


def test_fitting_brute():

    surf = read_surf(SMOOTH_FILE)
    pial = read_surf(PIAL_FILE, normals=False)

    n_rows = 8
    n_cols = 4
    rotation = 0
    start_label = 'elec008'
    start_vert = 31859

    grid = construct_grid(surf, start_vert, start_label, n_rows, n_cols, rotation=rotation)
    gamma = compute_distance(grid, pial)

    res = fitting_brute(surf, start_vert, start_label, n_rows, n_cols, gamma, pial)
    assert_array_almost_equal(
        res[0],
        array([0, 0, 0]),
        decimal=3)

    res = fitting_hop(surf, start_vert, start_label, n_rows, n_cols, gamma, pial)
    assert_array_almost_equal(
        res[0],
        array([0, 0, 0]),
        decimal=3)
