from numpy import array
from numpy.testing import assert_array_almost_equal

from gridgen.construct import construct_grid, make_grid
from gridgen.fitting import fitting_brute, fitting_hopping
from gridgen.io import read_surf
from gridgen.morphology.distance import compute_distance

from .paths import SMOOTH_FILE, PIAL_FILE


def notest_fitting_brute():

    surf = read_surf(SMOOTH_FILE)
    pial = read_surf(PIAL_FILE, normals=False)

    n_rows = 8
    n_cols = 4
    rotation = 0
    start_label = 'elec008'
    start_vert = 31859

    grid2d = make_grid(n_rows, n_cols)
    grid = construct_grid(surf, start_vert, start_label, grid2d, rotation=rotation)
    gamma = compute_distance(grid, pial)

    res = fitting_brute(surf, start_vert, start_label, grid2d, gamma, pial)
    assert_array_almost_equal(
        res[0],
        array([0, 0, 0]),
        decimal=3)


def notest_fitting_hop():

    surf = read_surf(SMOOTH_FILE)
    pial = read_surf(PIAL_FILE, normals=False)

    n_rows = 8
    n_cols = 4
    rotation = 0
    start_label = 'elec008'
    start_vert = 31859

    grid2d = make_grid(n_rows, n_cols)
    grid = construct_grid(surf, start_vert, start_label, grid2d, rotation=rotation)
    gamma = compute_distance(grid, pial)

    fitting_hopping(surf, start_vert, start_label, grid2d, gamma, pial)
