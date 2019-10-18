from numpy import array, cross
from numpy.linalg import norm


def calc_plane_to_axis(v):
    """Order of cross-product is very important.
    The normal of the vector is pointing towards the viewer, the first axis of
    the plane points up, so by the right-hand rule, the second axis needs to
    point right (so that it's comparable to the orientation of an array.)
    """
    perpendicular_v = array([0, 0, 1])
    v1 = cross(perpendicular_v, v)
    v1 /= norm(v1)  # rows
    v2 = cross(v, v1)  # columns
    v2 /= norm(v2)  # rows

    return array([v2, v1])
