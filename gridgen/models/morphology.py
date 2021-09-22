from numpy import dot, arccos, pi, zeros, cross, dtype, nanmin, NaN, errstate, ones, einsum
from numpy.linalg import norm
from scipy.stats import norm as normal_dist

from ..io import WIRE

d_ = dtype([
    ('label', '<U256'),   # labels cannot be longer than 256 char
    ('value', 'f4'),
    ])

EPSILON = 1e-5


def compute_morphology(grid, pial, distance='minimum', maximum_distance=None, penalty=1):
    if distance == 'ray':
        dist = _distance_ray(grid, pial)

    elif distance == 'minimum':
        dist = _distance_minimum(grid, pial)

    elif distance == 'view':
        dist = _distance_view(grid, pial)

    elif distance == 'cylinder':
        dist = _distance_cylinder(grid, pial)

    elif distance == 'pdf':
        dist = _distance_pdf(grid, pial)

    if maximum_distance and distance in ('ray', 'minimum'):
        with errstate(invalid='ignore'):
            i = dist['value'] > maximum_distance
        dist['value'][i] = NaN

    dist['value'] = dist['value'] ** (-1 * penalty)

    return dist


def _distance_ray(grid, pial):
    distance = zeros((grid.shape[0], grid.shape[1]), dtype=d_)
    distance['value'].fill(NaN)

    for i_x in range(grid.shape[0]):
        for i_y in range(grid.shape[1]):
            if grid['label'][i_x, i_y] != WIRE:
                dist = intersect_ray_triangle(
                    pial['pos'][pial['tri'][:, 0], :],
                    pial['pos'][pial['tri'][:, 1], :],
                    pial['pos'][pial['tri'][:, 2], :],
                    grid['pos'][i_x, i_y],
                    -1 * grid['norm'][i_x, i_y],
                    line=True
                    )

                distance['value'][i_x, i_y] = nanmin(dist)

    distance['label'] = grid['label']
    return distance


def _distance_minimum(grid, pial):
    distance = zeros((grid.shape[0], grid.shape[1]), dtype=d_)

    for i_x in range(grid.shape[0]):
        for i_y in range(grid.shape[1]):
            if grid['label'][i_x, i_y] == WIRE:
                continue
            distance['value'][i_x, i_y] = norm(pial['pos'] - grid['pos'][i_x, i_y], axis=1).min()

    distance['label'] = grid['label']
    return distance


def _distance_view(grid, pial):
    max_dist = 8
    max_angle = 15
    distance = zeros((grid.shape[0], grid.shape[1]), dtype=d_)
    distance['value'][:, :] = max_dist

    for i_x in range(grid.shape[0]):
        for i_y in range(grid.shape[1]):

            pos = grid['pos'][i_x, i_y]
            norm0 = grid['norm'][i_x, i_y]

            points = pial['pos'][norm(pial['pos'] - pos, axis=1) < max_dist, :]
            directions = (points - pos) / norm(points - pos, axis=1)[:, None]
            x = arccos(dot(directions, norm0 * -1)) / pi * 180
            if x.min() >= max_angle:
                continue
            distance['value'][i_x, i_y] = norm(points[x <= max_angle, :] - pos, axis=1).min()

    distance['label'] = grid['label']
    return distance


def _distance_cylinder(grid, pial):
    max_dist_to_elec = 100
    max_dist_to_line = 5
    distance = zeros((grid.shape[0], grid.shape[1]), dtype=d_)
    distance['value'][:, :] = max_dist_to_elec

    for i_x in range(grid.shape[0]):
        for i_y in range(grid.shape[1]):

            pos = grid['pos'][i_x, i_y]
            norm0 = grid['norm'][i_x, i_y]

            d = norm(pial['pos'] - pos, axis=1)
            # points = pial['pos'][d < max_dist_to_elec, :]
            dist_to_line = norm(cross(norm0, pos - pial['pos']), axis=1)
            distance['value'][i_x, i_y] = d[dist_to_line < max_dist_to_line].min()

    distance['label'] = grid['label']
    return distance


def _distance_pdf(grid, pial):
    distance = zeros((grid.shape[0], grid.shape[1]), dtype=d_)

    for i_x in range(grid.shape[0]):
        for i_y in range(grid.shape[1]):
            distance['value'][i_x, i_y] = normal_dist.pdf(norm(pial['pos'] - grid['pos'][i_x, i_y], axis=1), scale=2).sum()

    distance['label'] = grid['label']
    return distance


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
    i_good = ones(vertex0.shape[0], dtype=bool)
    edge1 = vertex1 - vertex0
    edge2 = vertex2 - vertex0
    h = cross(rayVector, edge2, axis=-1)
    a = einsum('ij,ij->i', edge1, h)  # dot product

    i_parallel = (a > -EPSILON) & (a < EPSILON)
    i_good[i_parallel] = False

    f = 1 / a

    s = rayOrigin - vertex0
    u = f * einsum('ij,ij->i', s, h)
    i_outside = (u < 0.0) | (u > 1.0)
    i_good[i_outside] = False

    q = cross(s, edge1)
    v = f * einsum('j,ij->i', rayVector, q)

    i_outside = (v < 0.0) | (u + v > 1.0)
    i_good[i_outside] = False

    t = f * einsum('ij,ij->i', edge2, q)
    if not line:
        i_opposite = t < EPSILON
        i_good[i_opposite] = NaN

    t[~i_good] = NaN

    return t
