from logging import getLogger
from numpy import array, arange, argmin, cross, pi, sqrt
from numpy.linalg import norm
from scipy.spatial.transform import Rotation

from .algebra import calc_plane_to_axis

interelec_distance = 3
MAX_ANGLE = 10
POSSIBLE_DEGREES = arange(-MAX_ANGLE, MAX_ANGLE)

lg = getLogger(__name__)


def find_new_pos_0d(x, y, grid, neighbors, surf, radians=0):

    x_neighbor, y_neighbor = neighbors[x, y, 0, :]
    pos_neighbor = grid[x_neighbor, y_neighbor, 0, :]  # this cannot be nan
    normal_neighbor = grid[x_neighbor, y_neighbor, 1, :]  # this cannot be nan

    coords_2d = array([x - x_neighbor, y - y_neighbor]) * interelec_distance

    plane = calc_plane_to_axis(normal_neighbor, radians=radians)
    if coords_2d[0] == 0:
        rotation_axis = plane[0, :]
    else:
        rotation_axis = plane[1, :]

    pos_potential = []
    plane_potential = []
    distance = []
    for degrees in POSSIBLE_DEGREES:

        r = Rotation.from_rotvec(rotation_axis * degrees / 180 * pi)
        trans_2d_to_3d = plane @ r.as_dcm()
        new_pos = coords_2d @ trans_2d_to_3d + pos_neighbor
        pos_potential.append(new_pos)
        plane_potential.append(trans_2d_to_3d)
        distance.append(norm(surf['pos'] - new_pos, axis=1).min())

    idx_min_angle = argmin(distance)
    min_angle = POSSIBLE_DEGREES[idx_min_angle]

    new_pos = pos_potential[idx_min_angle]
    new_plane = plane_potential[idx_min_angle]
    new_normal = cross(new_plane[0, :], new_plane[1, :])
    lg.debug(f'New point in grid row: {x}, column: {y}')
    lg.debug(f'\tpos: {new_pos}\n\tnormal: {new_normal}')
    lg.info(f'Minimum angle {min_angle}°, distance to surface {min(distance):.3f}')

    return new_pos, new_normal


def find_new_pos_2d(x, y, grid, neighbors, surf):
    """Make line from first neighbor to second neighbor, and the new point
    should be on the right of that line
    """
    x1, y1 = neighbors[x, y, 0, :]
    x2, y2 = neighbors[x, y, 1, :]

    pos1 = grid[x1, y1, 0, :]
    pos2 = grid[x2, y2, 0, :]
    normal1 = grid[x1, y1, 1, :]
    normal2 = grid[x2, y2, 1, :]

    rotation_axis = pos1 - pos2
    rotation_axis /= norm(rotation_axis)

    center = (pos1 + pos2) / 2
    normal_center = (normal1 + normal2) / 2
    normal_center /= norm(normal_center)  # avoid rounding errors

    search_direction = cross(normal_center, rotation_axis)

    # this will fail when the distance between two opposite electrodes is more than twice the interelectrode distance
    search_distance = sqrt((interelec_distance ** 2 - norm(center - pos1) ** 2))

    pos_potential = []
    distance = []
    for degrees in POSSIBLE_DEGREES:

        r = Rotation.from_rotvec(rotation_axis * degrees / 180 * pi)
        new_pos = (search_distance * search_direction) @ r.as_dcm() + center
        pos_potential.append(new_pos)
        distance.append(norm(surf['pos'] - new_pos, axis=1).min())

    idx_min_angle = argmin(distance)
    min_angle = POSSIBLE_DEGREES[idx_min_angle]
    new_pos = pos_potential[idx_min_angle]
    r = Rotation.from_rotvec(rotation_axis * min_angle / 180 * pi)
    new_normal = normal_center @ r.as_dcm()

    lg.debug(f'New point in grid row: {x}, column: {y}')
    lg.debug(f'\tpos: {new_pos}\n\tnormal: {new_normal}')
    lg.info(f'Minimum angle {min_angle}°, distance to surface {min(distance):.3f}')

    return new_pos, new_normal
