from logging import getLogger
from numpy import array, arange, argmin, cross, pi
from numpy.linalg import norm
from scipy.spatial.transform import Rotation

from .algebra import calc_plane_to_axis

interelec_distance = 3

lg = getLogger(__name__)


def find_new_pos_1d(x, y, grid, neighbors, surf):

    x_neighbor, y_neighbor = neighbors[x, y, 0, :]
    pos_neighbor = grid[x_neighbor, y_neighbor, 0, :]  # this cannot be nan
    normal_neighbor = grid[x_neighbor, y_neighbor, 1, :]  # this cannot be nan

    coords_2d = array([x - x_neighbor, y - y_neighbor]) * interelec_distance

    plane = calc_plane_to_axis(normal_neighbor)
    if coords_2d[0] == 0:
        rotation_axis = plane[0, :]
    else:
        rotation_axis = plane[1, :]

    pos_potential = []
    plane_potential = []
    distance = []
    possible_degrees = arange(-30, 30)
    for degrees in possible_degrees:

        r = Rotation.from_rotvec(rotation_axis * degrees / 180 * pi)
        trans_2d_to_3d = plane @ r.as_dcm()
        new_pos = coords_2d @ trans_2d_to_3d + pos_neighbor
        pos_potential.append(new_pos)
        plane_potential.append(trans_2d_to_3d)
        distance.append(norm(surf['pos'] - new_pos, axis=1).min())

    idx_min_angle = argmin(distance)
    min_angle = possible_degrees[idx_min_angle]

    new_pos = pos_potential[idx_min_angle]
    new_plane = plane_potential[idx_min_angle]
    new_normal = cross(new_plane[0, :], new_plane[1, :])
    lg.debug(f'New point in grid row: {x}, column: {y}')
    lg.debug(f'\tpos: {new_pos}\n\tnormal: {new_normal}')
    lg.info(f'Minimum angle {min_angle}, distance to surface {min(distance)}')

    return new_pos, new_normal
