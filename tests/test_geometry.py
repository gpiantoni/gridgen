from numpy import pi, array, dot
from numpy.testing import assert_array_almost_equal
from gridloc.geometry import calc_plane_to_axis


def test_rotation():
    v = array([5., 1., 2.])
    v1 = calc_plane_to_axis(v, radians=0)
    v2 = calc_plane_to_axis(v, radians=pi / 2)

    assert_array_almost_equal(
        v1[0, :],
        v2[1, :])
    assert dot(v1[0, :], v2[0, :]) < 1e-6
    assert dot(v1[1, :], v2[1, :]) < 1e-6
