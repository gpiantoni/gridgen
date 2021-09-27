from numpy.testing import assert_almost_equal
from .paths import T1_FILE, SMOOTH_FILE, PIAL_FILE, ECOG_FILE
from .paths import EXAMPLES as parameters

from gridgen.bin.parameters import validate_template, TEMPLATE
from gridgen.ecog import read_ecog, put_ecog_on_grid2d
from gridgen.models import compare_model_with_ecog, sum_models, make_grid3d_model
from gridgen.models.morphology import compute_morphology
from gridgen.io import read_surface_ras_shift, read_surf, read_mri
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


def test_combine():

    grid2d = make_grid_with_labels(4, 5, 'TBLR', 'chan{}')

    mris = read_mri(**parameters['mri'])
    for k in ('grid3d', 'functional', 'morphology', 'fit'):
        parameters[k] = validate_template(TEMPLATE[k], parameters.get(k, {}))

    model = make_grid3d_model(mris, grid2d, parameters['grid3d'], parameters['initial'], parameters['morphology'], parameters['functional'])

    parameters['fit']['metric'] = 'sum'

    m1 = model.copy()
    m1['functional'] = None
    out = sum_models(m1, parameters['fit'])[1]
    assert_almost_equal(out, 7.434, decimal=3)

    m2 = model.copy()
    m2['morphology'] = None
    out = sum_models(m2, parameters['fit'])[1]
    assert_almost_equal(out, 59438.953, decimal=3)

    del parameters['fit']['functional_contribution']
    out = sum_models(model, parameters['fit'])[1]
    assert_almost_equal(out, 53306.891, decimal=3)

    tf = read_ecog(ECOG_FILE)
    ecog2d = put_ecog_on_grid2d(tf, grid2d)

    parameters['fit']['metric'] = 'nonparametric'
    out = compare_model_with_ecog(m1, ecog2d, parameters['fit'])[1]
    assert_almost_equal(out, 0.264, decimal=3)

    compare_model_with_ecog(model, ecog2d, parameters['fit'])
