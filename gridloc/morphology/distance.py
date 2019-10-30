from numpy import dot, arccos, pi, ones, zeros, norm
from scipy.stats import norm as normal_dist


def compute_distance(grid, pial, method='minimum'):

    if method == 'minimum':
        distance = _distance_minimum(grid, pial)

    elif method == 'view':
        _distance_view(grid, pial)

    elif method == 'pdf':
        _distance_view(grid, pial)

    return distance


def _distance_minimum(grid, pial):
    distance = zeros((grid.shape[0], grid.shape[1]))

    for i_x in range(grid.shape[0]):
        for i_y in range(grid.shape[1]):
            distance[i_x, i_y] = norm(pial['pos'] - grid[i_x, i_y, 0, :], axis=1).min()

    return distance


def _distance_view(grid, pial):
    max_dist = 8
    max_angle = 15
    distance = ones((grid.shape[0], grid.shape[1])) * max_dist

    for i_x in range(grid.shape[0]):
        for i_y in range(grid.shape[1]):

            pos = grid[i_x, i_y, 0, :]
            norm0 = grid[i_x, i_y, 1, :]

            points = pial['pos'][norm(pial['pos'] - pos, axis=1) < max_dist, :]
            directions = (points - pos) / norm(points - pos, axis=1)[:, None]
            x = arccos(dot(directions, norm0 * -1)) / pi * 180
            if x.min() >= max_angle:
                continue
            distance[i_x, i_y] = norm(points[x <= max_angle, :] - pos, axis=1).min()

    return distance


def _distance_pdf(grid, pial):
    distance = zeros((grid.shape[0], grid.shape[1]))

    for i_x in range(grid.shape[0]):
        for i_y in range(grid.shape[1]):
            distance[i_x, i_y] = normal_dist.pdf(norm(pial['pos'] - grid[i_x, i_y, 0, :], axis=1), scale=2).sum()

    return distance
