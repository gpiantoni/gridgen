from numpy import cross
from numpy.linalg import norm
from numpy.testing import assert_array_almost_equal

from gridloc.matlab.geometry import project_to_cortex
from gridloc.io import read_surf

from .paths import PIAL_FILE


def test_project_raytriangle():
    surf = read_surf(PIAL_FILE, normals=False)
    i_pos = 9099
    distance = 2

    # reconstruct normal for one plane (it's not necessary to compute them all)
    tris = surf['pos'][surf['tri'][i_pos, :]]
    tri_norm = cross(tris[:, 1] - tris[:, 0], tris[:, 2] - tris[:, 0])
    tri_norm /= norm(tri_norm)

    # create fake point
    point = tris.mean(axis=0) + tri_norm * distance
    direction = tri_norm * -1

    out_distance, out_point = project_to_cortex(surf, point, direction)

    assert_array_almost_equal(distance, out_distance)
    assert_array_almost_equal(tris.mean(axis=0), out_point)
