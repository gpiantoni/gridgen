"""Functions to compute vector-based transformation and searches.
"""
from numpy import array, zeros, argmin, cross
from numpy.linalg import norm
from scipy.spatial.transform import Rotation


NEIGHBORS = array([
    [-1, 0],
    [0, 1],
    [1, 0],
    [0, -1],
    ])


def count_neighbors(grid, x, y):
    """Count neighbors around one electrode

    Parameters
    ----------
    grid : 2d ndarray
        grid (n_rows, n_columns) with fields (label, pos, norm, done)
    x : int
        row index for this electrode
    y : int
        column index for this electrode

    Returns
    -------
    int
        number of neighbors
    list of tuple
        index of the direct neighbor to this channel
    """
    n_rows, n_cols = grid.shape

    idx = zeros(4, bool)

    i = 0
    for n_x, n_y in NEIGHBORS:
        if (0 <= (x + n_x) < n_rows) and (0 <= (y + n_y) < n_cols) and grid['done'][x + n_x, y + n_y]:
            idx[i] = True
        i += 1

    n = idx.sum()

    """make sure that any indexing does not have more than 2 neighbors
    Also, when you have two neighbors, the order is important. If you draw a
    line from the first neighbor to the second neighbor looking from above,
    the target position should be on the right"""
    if (idx == array([True, False, False, True])).all():
        coords = NEIGHBORS[(3, 0), :] + (x, y)
    else:
        coords = NEIGHBORS[idx] + (x, y)
    return n, coords


def search_grid(surf, ref_vert, x, y):
    """Given a reference vertex, compute a plane perpendicular to the normal of
    that point. Then, move x-mm in the x direction and y-mm in the y direction,
    then find the closest vertex to that point.

    Parameters
    ----------
    surf : dict
        surface with normals. It's better to use the smooth surface to have
        reasonable results
    ref_vert : int
        index of the reference vertex (the plane will be computed relative to
        the normal of this vertex)
    x : float
        distance in mm from the reference vertex in one direction
    y : float
        distance in mm from the reference vertex in the direction perpendicular
        to x

    Returns
    -------
    int
        index of the vertex which is roughly x-mm and y-mm away

    Notes
    -----
    Because vertices are not regularly distributed on the surface, the distance
    between the reference vertex and the output vertex does not need to be equal
    to `sqrt(x ** 2 + y ** 2)`, but it should be close enough.

    There are an infinite number of planes perpendicular to the normal. x-y
    define one of the possible planes. The plane is always the same for one
    normal (due to internal convention). You can rotate the plane later.
    """
    pos = surf['pos'][ref_vert, :]
    normal = surf['pos_norm'][ref_vert, :]

    coords_2d = array([x, y])
    plane = calc_plane_to_axis(normal)
    target = coords_2d @ plane + pos
    new_vector = argmin(norm(surf['pos'] - target, axis=1))

    return new_vector


def calc_plane_to_axis(v, radians=0):
    """Compute the plane perpendicular to the input vector. Because there are
    infinite number of vectors belonging to one plane, we fix the first vector
    to point towards the superior part of the brain.

    Parameters
    ----------
    v : array of 3 values
        input vector. The dot product of this vector and each of the two output
        vectors is 0
    radians : float
        rotation to apply about the input v vector, clockwise.

    Returns
    -------
    2x3 array
        2 vectors defining the plane perpendicular to input vector

    Notes
    -----
    Order of cross-product is very important.
    The normal of the vector is pointing towards the viewer, the first axis of
    the plane points up, so by the right-hand rule, the second axis needs to
    point right (so that it's comparable to the orientation of an array.)
    """
    v /= norm(v)  # essential, otherwise it's not correct
    perpendicular_v = array([0., 0., 1.])
    v1 = cross(perpendicular_v, v)
    v1 /= norm(v1)  # rows
    v2 = cross(v1, v)  # columns
    v2 /= norm(v2)  # rows

    plane = array([v2, v1])

    # apply rotation about normal (normal points to you, then it's clockwise)
    r = Rotation.from_rotvec(v * radians)
    return plane @ r.as_matrix()
