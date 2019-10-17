from numpy import array, abs, cross
from numpy.linalg import norm


def calc_plane_to_axis(v):
    # TODO: choose axis 1 and 0 based on the main direction of v
    v1 = array([1, 0, abs(v[0] / v[2])])  # make sure it's always up
    v1 = v1 / norm(v1)
    v2 = cross(v, v1)
    return array([v1, v2])
