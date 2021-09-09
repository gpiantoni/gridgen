"""Sub-functions to find one node (electrode) at the time, based on the location
of the previous nodes"""
from logging import getLogger
from numpy import array, arange, argmin, cross, pi, sqrt
from numpy.linalg import norm
from scipy.spatial.transform import Rotation

from .geometry import calc_plane_to_axis

# POSSIBLE_DEGREES = arange(-MAX_ANGLE, MAX_ANGLE, 0.2)

lg = getLogger(__name__)
lg.propagate = False


def find_new_pos_0d(grid, neighbors, surf, x, y, radians=0):
    """Find the position of an electrode which has 1 neighbor perpendicular to it.
    We try to estimate the position of electrode `x` based on the position of `a`.

        a -- x

    In this example, radians would be π/2, because `x` is at the right side of `a`.

    Parameters
    ----------
    grid : 2d ndarray
        grid (n_rows, n_columns) with fields (label, pos, norm, done)
    neighbors : 2d array
        the row specifies the x and y coordinates of the neighbor `a`.
    surf : dict
        surface with normals. It's better to use the smooth surface to have
        reasonable results
    x : int
        row index for this electrode
    y : int
        column index for this electrode
    radians : float
        rotation from the inferior-superior axis (clockwise, in radians)

    Returns
    -------
    2d ndarray
        grid (n_rows, n_columns) with fields (label, pos, norm, done), with the
        position and normal of the new electrode included
    """
    x_neighbor, y_neighbor = neighbors[0, :]
    pos_neighbor = grid['pos'][x_neighbor, y_neighbor]  # this cannot be nan
    normal_neighbor = grid['norm'][x_neighbor, y_neighbor]  # this cannot be nan

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
        trans_2d_to_3d = plane @ r.as_matrix()
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

    grid['pos'][x, y] = new_pos
    grid['norm'][x, y] = new_normal
    grid['done'][x, y] = True


def find_new_pos_1d(grid, neighbors, surf, x, y, opposite):
    """Find the position of an electrode which has 1 neighbor on the same line.
    We try to estimate the position of electrode `x` based on the position of `a`
    and `b`.

        a -- b -- x

    Parameters
    ----------
    grid : 2d ndarray
        grid (n_rows, n_columns) with fields (label, pos, norm, done)
    neighbors : 2d array
        the row specifies the x and y coordinates of the neighbor `b`.
    surf : dict
        surface with normals. It's better to use the smooth surface to have
        reasonable results
    x : int
        row index for this electrode
    y : int
        column index for this electrode
    opposite : tuple of 2 floats
        x and y indices of the electrode opposite to the electrode of interest
        (i.e. `a`)

    Returns
    -------
    2d ndarray
        grid (n_rows, n_columns) with fields (label, pos, norm, done), with the
        position and normal of the new electrode included
    """
    x1, y1 = neighbors[0, :]
    x2, y2 = opposite

    pos1 = grid['pos'][x1, y1]
    pos2 = grid['pos'][x2, y2]
    normal1 = grid['norm'][x1, y1]

    direction = (pos1 - pos2)
    direction /= norm(direction)
    rotation_axis = cross(direction, normal1)

    pos_potential = []
    plane_potential = []
    distance = []
    for degrees in POSSIBLE_DEGREES:

        r = Rotation.from_rotvec(rotation_axis * degrees / 180 * pi)
        new_pos = direction @ r.as_matrix() * interelec_distance + pos1
        new_normal = cross(rotation_axis, direction @ r.as_matrix())
        pos_potential.append(new_pos)
        plane_potential.append(new_normal)
        distance.append(norm(surf['pos'] - new_pos, axis=1).min())

    idx_min_angle = argmin(distance)
    min_angle = POSSIBLE_DEGREES[idx_min_angle]

    new_pos = pos_potential[idx_min_angle]
    r = Rotation.from_rotvec(rotation_axis * min_angle / 180 * pi)
    new_normal = cross(rotation_axis, direction @ r.as_matrix())

    lg.debug(f'New point in grid row: {x}, column: {y}')
    lg.debug(f'\tpos: {new_pos}\n\tnormal: {new_normal}')
    lg.info(f'Minimum angle {min_angle}°, distance to surface {min(distance):.3f}')

    grid['pos'][x, y] = new_pos
    grid['norm'][x, y] = new_normal
    grid['done'][x, y] = True


def find_new_pos_2d(grid, neighbors, surf, x, y):
    """Find the position of an electrode which has already 2 neighbors. So, we
    try to estimate the position of electrode `x` based on the position of `a`
    and `b`.

        b -- x
        |    |
        |    |
        c -- a

    Make a line from `a` to `b`, and the new point `x` should be on the right
    of that line (convention). Note that point `c` is not included in calculations

    Parameters
    ----------
    grid : 2d ndarray
        grid (n_rows, n_columns) with fields (label, pos, norm, done)
    neighbors : 2d array
        each row specifies the x and y coordinates of the two neighbors. The
        order is very important in this case and it should have been taken
        care of by `count_neighbors()`.
    surf : dict
        surface with normals. It's better to use the smooth surface to have
        reasonable results
    x : int
        row index for this electrode
    y : int
        column index for this electrode

    Returns
    -------
    2d ndarray
        grid (n_rows, n_columns) with fields (label, pos, norm, done), with the
        position and normal of the new electrode included
    """
    x1, y1 = neighbors[0, :]
    x2, y2 = neighbors[1, :]

    pos1 = grid['pos'][x1, y1]
    pos2 = grid['pos'][x2, y2]
    normal1 = grid['norm'][x1, y1]
    normal2 = grid['norm'][x2, y2]

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
        new_pos = (search_distance * search_direction) @ r.as_matrix() + center
        pos_potential.append(new_pos)
        distance.append(norm(surf['pos'] - new_pos, axis=1).min())

    idx_min_angle = argmin(distance)
    min_angle = POSSIBLE_DEGREES[idx_min_angle]
    new_pos = pos_potential[idx_min_angle]
    r = Rotation.from_rotvec(rotation_axis * min_angle / 180 * pi)
    new_normal = normal_center @ r.as_matrix()

    lg.debug(f'New point in grid row: {x}, column: {y}')
    lg.debug(f'\tpos: {new_pos}\n\tnormal: {new_normal}')
    lg.info(f'Minimum angle {min_angle}°, distance to surface {min(distance):.3f}')

    grid['pos'][x, y] = new_pos
    grid['norm'][x, y] = new_normal
    grid['done'][x, y] = True
