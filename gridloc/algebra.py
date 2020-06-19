from numpy import array, cross
from numpy.linalg import norm

from scipy.spatial.transform import Rotation


def calc_plane_to_axis(v, radians=0):
    """Compute the plane perpendicular to the input vector. Because there are
    infinite number of vectors belonging to one plane, we fix the first vector
    to point towards the superior part of the brain.

    Parameters
    ----------
    v : array of 3 values
        input vector. The dot product of this vector and each of the two output
        vectors is 0
    radians : float
        rotation to apply about the input v vector, clockwise.

    Returns
    -------
    2x3 array
        2 vectors defining the plane perpendicular to input vector

    Notes
    -----
    Order of cross-product is very important.
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
    r = Rotation.from_rotvec(v * radians)
    return plane @ r.as_matrix()
