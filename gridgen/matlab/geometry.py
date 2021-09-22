from numpy import NaN, isnan, nanargmin, array, dot, sum, where
from numpy.linalg import norm

from ..models.morphology import intersect_ray_triangle


def project_to_cortex(surf, point, normal, sorted_triangles=None):
    """Project a point (electrode) onto the triangulated mesh (surface).

    Parameters
    ----------
    surf : dict with 'pos', 'tri'
        surface of the brain use to project the electrodes (it's not necessary
        to have 'tri_norm')
    point : (3, ) array
        x, y, z coordinates of the electrodes
    normal : (3, ) array
        normal of the electrodes
    sorted_triangles : (n, ) array
        indices of triangles to select (order is important). If you pass this
        parameters, it will return the first triangle which is intersected by
        the line

    Returns
    -------
    int
        distance between surface and electrode
    array (3, )
        projected position onto the surface


    Notes
    -----
    Returns NaN values when there is no intersection possible to the cortex.
    """
    if sorted_triangles is None:
        vertices = surf['tri']
    else:
        vertices = surf['tri'][sorted_triangles]

    normal = normal / norm(normal)
    t = intersect_ray_triangle(
        surf['pos'][vertices][:, 0, :],
        surf['pos'][vertices][:, 1, :],
        surf['pos'][vertices][:, 2, :],
        point,
        normal,
        line=True)

    if isnan(t).all():
        return NaN, array([NaN, NaN, NaN])

    if sorted_triangles is None:
        i = nanargmin(abs(t))
    else:
        i = where(~isnan(t))[0][0]
    projected_point = point + normal * t[i]

    return t[i], projected_point


def check_if_point_in_triangle(v0, v1, v2, p):
    """They have to be on the same plane"""
    p_bary = cartesian_to_barycentric(v0, v1, v2, p)
    return (0 <= p_bary).all() and (p_bary <= 1).all()


def intersect_line_plane(vertex0, tri_norm, pos, pos_norm):
    pos_norm = pos_norm / norm(pos_norm)
    t = (dot(tri_norm, vertex0) - dot(tri_norm, pos)) / dot(pos_norm, tri_norm)
    intersection = pos + pos_norm * t
    return t, intersection


def cartesian_to_barycentric(v0, v1, v2, p):
    S = sum((v0 - v1) * (v0 - v2))
    S1 = sum((v1 - p) * (v2 - p))
    S2 = sum((p - v2) * (p - v1))

    a = S1 / S
    b = S2 / S
    c = 1 - a - b
    return array([a, b, c])
