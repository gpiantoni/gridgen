from numpy import array, cross
from numpy.linalg import norm

from scipy.spatial.transform import Rotation


def calc_plane_to_axis(v, degrees=0):
    """Order of cross-product is very important.
    The normal of the vector is pointing towards the viewer, the first axis of
    the plane points up, so by the right-hand rule, the second axis needs to
    point right (so that it's comparable to the orientation of an array.)
    """
    v /= norm(v)  # essential, otherwise it's not correct
    perpendicular_v = array([0., 0., 1.])
    v1 = cross(perpendicular_v, v)
    v1 /= norm(v1)  # rows
    v2 = cross(v1, v)  # columns
    v2 /= norm(v2)  # rows

    plane = array([v2, v1])

    # apply rotation about normal (normal points to you, then it's clockwise)
    r = Rotation.from_rotvec(v * degrees)
    return plane @ r.as_dcm()
