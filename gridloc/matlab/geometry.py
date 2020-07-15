from numpy import cross, NaN, einsum, empty, isnan, nanargmin
from numpy.linalg import norm

EPSILON = 1e-5


def project_to_cortex(surf, point, normal):
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

    Returns
    -------
    int
        distance between surface and electrode
    array (3, )
        projected position onto the surface
    """
    normal /= norm(normal)
    t = intersect_ray_triangle(
        surf['pos'][surf['tri'][:, 0]],
        surf['pos'][surf['tri'][:, 1]],
        surf['pos'][surf['tri'][:, 2]],
        point,
        normal)

    i = nanargmin(t)
    projected_point = point + normal * t[i]

    return t[i], projected_point


def intersect_ray_triangle(vertex0, vertex1, vertex2, rayOrigin, rayVector):
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
    i_opposite = t < EPSILON
    out_dist[i_opposite] = NaN

    i_good = ~isnan(out_dist)
    out_dist[i_good] = t[i_good]

    return out_dist
