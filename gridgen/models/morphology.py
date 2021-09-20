from numpy import dot, arccos, pi, zeros, cross, dtype, nanmin, NaN, isnan
from numpy.linalg import norm
from scipy.stats import norm as normal_dist

from ..matlab.geometry import intersect_ray_triangle
from ..io import WIRE

d_ = dtype([
    ('label', '<U256'),   # labels cannot be longer than 256 char
    ('value', 'f4'),
    ])


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
        i = isnan(dist['value']) | (dist['value'] > maximum_distance)
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
            points = pial['pos'][d < max_dist_to_elec, :]
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
