from numpy import array, arange, argmin, cross, pi
from numpy.linalg import norm
from scipy.spatial.transform import Rotation

from .algebra import calc_plane_to_axis

interelec_distance = 3

def _find_new_pos_1d(x, y, grid, neighbors, surf):

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
    possible_degrees = arange(-30, 30)
    for degrees in possible_degrees:

        r = Rotation.from_rotvec(rotation_axis * degrees / 180 * pi)
        pos_potential.append(
            coords_2d @ plane @ r.as_dcm() + pos_neighbor
            )
        plane_potential.append(
            plane @ r.as_dcm()
            )

    idx_min_angle = argmin(array([norm(surf['pos'] - x, axis=1).min() for x in pos_potential]))
    min_angle = possible_degrees[idx_min_angle]

    print(min_angle)

    new_pos = pos_potential[idx_min_angle]
    new_plane = plane_potential[idx_min_angle]

    new_normal = cross(new_plane[0, :], new_plane[1, :])
    return new_pos, new_normal
