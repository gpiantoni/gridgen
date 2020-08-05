from numpy import cross, NaN, einsum, empty, isnan, nanargmin, array, dot, sum, where
from numpy.linalg import norm

EPSILON = 1e-5


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
        i = nanargmin(t)
    else:
        i = where(~isnan(t))[0]
    projected_point = point + normal * t[i]

    return t[i], projected_point


def intersect_ray_triangle(vertex0, vertex1, vertex2, rayOrigin, rayVector, line=False):
    """Implementation of the Möller-Trumbore algorithm to detect ray-triangle
    intersection.


    Parameters
    ----------
    vertex0 : (n, 3) array
        the first vertex for all the triangles of the mesh
    vertex1 : (n, 3) array
        the second vertex for all the triangles of the mesh
    vertex2 : (n, 3) array
        the third vertex for all the triangles of the mesh
    rayOrigin : (3, ) array
        start position of the ray
    rayVector : (3, ) array
        direction of the ray vector
    line : bool
        if False, the sign of the normal is important (treated as a ray, leaving rayOrigin).
        if True, it detects intersection in both direction (treated as a line
        with two directions)

    Returns
    -------
    array (n, )
        distance of each triangle to the point of interest. It's NaN if 1.
        the direction of the point is parallel to the triangle, 2. the projected
        point is outside the triangle, 3. the direction is in the opposite direction
        to the triangle (so, this is a ray, not a line).

    Notes
    -----
    Implementation derived from https://en.wikipedia.org/wiki/M%C3%B6ller%E2%80%93Trumbore_intersection_algorithm
    """
    out_dist = empty(vertex0.shape[0])
    edge1 = vertex1 - vertex0
    edge2 = vertex2 - vertex0
    h = cross(rayVector, edge2, axis=-1)
    a = einsum('ij,ij->i', edge1, h)  # dot product

    i_parallel = (a > -EPSILON) & (a < EPSILON)
    out_dist[i_parallel] = NaN

    f = 1 / a

    s = rayOrigin - vertex0
    u = f * einsum('ij,ij->i', s, h)
    i_outside = (u < 0.0) | (u > 1.0)
    out_dist[i_outside] = NaN

    q = cross(s, edge1)
    v = f * einsum('j,ij->i', rayVector, q)

    i_outside = (v < 0.0) | (u + v > 1.0)
    out_dist[i_outside] = NaN

    t = f * einsum('ij,ij->i', edge2, q)
    if not line:
        i_opposite = t < EPSILON
        out_dist[i_opposite] = NaN

    i_good = ~isnan(out_dist)
    out_dist[i_good] = t[i_good]

    return out_dist


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
