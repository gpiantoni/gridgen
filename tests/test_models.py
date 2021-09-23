from numpy.testing import assert_almost_equal
from .paths import T1_FILE, SMOOTH_FILE, PIAL_FILE

from gridgen.models.morphology import compute_morphology
from gridgen.io import read_surface_ras_shift, read_surf
from gridgen.grid2d import make_grid_with_labels
from gridgen.grid3d import construct_grid


def test_morphology():

    offset = read_surface_ras_shift(T1_FILE)
    smooth = read_surf(SMOOTH_FILE, ras_shift=offset, normals=True)
    pial = read_surf(PIAL_FILE, ras_shift=offset, normals=False)

    grid2d = make_grid_with_labels(4, 3, 'TBLR', chan_pattern='elec{}')
    grid3d = {
        'interelec_distance': 3,
        'maximum_angle': 5,
        'step_angle': 0.2,
        }
    grid = construct_grid(smooth, 37613, 'elec1', grid2d['label'], grid3d, rotation=20)

    val = compute_morphology(grid, pial, distance='ray', maximum_distance=10, penalty=2)
    assert_almost_equal(
        val['value'][1, 0],
        0.071,
        decimal=3)

    val = compute_morphology(grid, pial, distance='minimum', maximum_distance=None, penalty=1)
    assert_almost_equal(
        val['value'][1, 0],
        0.407,
        decimal=3)

    val = compute_morphology(grid, pial, distance='cylinder', maximum_distance=5, penalty=1)
    assert_almost_equal(
        val['value'][1, 0],
        0.407,
        decimal=3)

    val = compute_morphology(grid, pial, distance='view', maximum_distance=5, penalty=1)
    assert_almost_equal(
        val['value'][1, 0],
        0.317,
        decimal=3)

    val = compute_morphology(grid, pial, distance='pdf', maximum_distance=5, penalty=2)
    assert_almost_equal(
        val['value'][1, 0],
        2.430,
        decimal=3)
