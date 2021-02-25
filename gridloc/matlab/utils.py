"""Functions directly useful to compat.py
"""

from os import nice
from numpy import arange, meshgrid, c_, zeros, prod, argmax, dot, cross, NaN, array, eye, ones, where, argsort, concatenate
from numpy.linalg import norm, solve
from scipy.spatial.transform import Rotation

EPSILON = 1e-5


def calcCoords(c, GridSteps, dims):
    """Calculate the coordinates of a grid placed with the point c as
    middlepoint of the space where the grid is going to be projected.
    Also define the GridSteps and GridSize for the density and
    size, respectively, of the mesh that's being calculated.
    """
    x_steps = arange(dims[0]) * GridSteps[0]
    y_steps = arange(dims[1]) * GridSteps[1]

    x_mesh, y_mesh = meshgrid(x_steps - x_steps.mean(), y_steps - y_steps.mean())
    diff_mat = c_[zeros(prod(dims)), x_mesh.flatten('F'), y_mesh.flatten('F')]

    return c + diff_mat


def plane_intersect(N1, A1, N2, A2):
    """Python implementation of plane_intersect.m

    Parameters
    ----------
    N1 : (3, ) array
        equation for the first plane
    A1 : (3, ) array
        point belonging to the first plane
    N2 : (3, ) array
        equation for the second plane
    A2 : (3, ) array
        point belonging to the second plane

    Returns
    -------
    array (3, )
        one point belonging to the intersection of the plane
    array (3, )
        equation describing the plane

    Notes
    -----
    It will return NaN if there is no intersection between the two planes.
    """
    P = array([0., 0., 0.])
    N = cross(N1, N2)
    if norm(N) < EPSILON:  # parallel or coincide
        return array([NaN, NaN, NaN]), NaN

    maxc = argmax(abs(N))

    d1 = -dot(N1, A1)
    d2 = -dot(N2, A2)

    if maxc == 0:
        P[0] = 0
        P[1] = (d2 * N1[2] - d1 * N2[2]) / N[0]
        P[2] = (d1 * N2[1] - d2 * N1[1]) / N[0]

    elif maxc == 1:
        P[0] = (d1 * N2[2] - d2 * N1[2]) / N[1]
        P[1] = 0
        P[2] = (d2 * N1[0] - d1 * N2[0]) / N[1]

    elif maxc == 2:
        P[0] = (d2 * N1[1] - d1 * N2[1]) / N[2]
        P[1] = (d1 * N2[0] - d2 * N1[0]) / N[2]
        P[2] = 0

    return P, N


def AxelRot(radians, u, x0):
    """Python implementation of AxelRot.m

    Parameters
    ----------
    radians : float
        rotation in radians, clockwise
    u : (3, ) array
        axis about to which compute the rotation
    x0 : (3, ) array
        point for shift

    Returns
    -------
    array (4, 4)
        affine matrix with rotation and translation
    """
    u = u / norm(u)
    AxisShift = x0 - (x0 @ u) * u  # l. 85

    Mshift = eye(4)
    Mshift[:3, 3] = -AxisShift

    Mroto = eye(4)
    Mroto[:3, :3] = Rotation.from_rotvec(u * radians).as_matrix()  # l.92

    # be careful about the sign of AxisShift
    return solve(Mshift, Mroto) @ Mshift  # l. 94


def _apply_affine(c, affine):
    """Apply affine matrix to many points at the same time.

    Parameters
    ----------
    c : (n, 3) array
        points to transform
    affine : (4, 4) array
        affine matrix (bottom row should be [0, 0, 0, 1])

    Returns
    -------
    array (n, 3)
        points transformed
    """
    C = c_[c, ones(c.shape[0])]
    X = (affine @ C.T).T
    return X[:, :3]


def _sort_closest_triangles(surf, electrode, intersval):
    """sort triangles based on the distance of the vertices. First compute
    which vertices are closer than `intersval` to the electrode, then select
    only vertices which have that vertex. This function keeps the order of distance
    of the vertices. Each vertex has three triangles, so the order of the first
    3 triangles is arbitrary. Here we use matlab convention for compatibility.

    Parameters
    ----------
    surf : dict with 'pos', 'tri'
        surface of the brain use to project the electrodes (it's not necessary
        to have 'tri_norm')
    electrode : (3, ) array
        x, y, z coordinates of the electrodes
    intersval : float
        radius in which the vertices should be

    Returns
    -------
    list
        sorted list of indices of triangles, sorted mostly by the distance of
        their closest vertex to `electrode`.
    """
    dvect = norm(electrode - surf['pos'], axis=1)
    closevert = where(dvect < intersval)[0]
    dvecti = argsort(dvect[closevert])
    sortvert = closevert[dvecti]

    # l. 176-192
    sorttri = []
    tri = surf['tri'].copy()
    for cv in sortvert:
        rows = concatenate([where(cv == tri[:, i])[0] for i in range(3)])
        tri[rows, :] = 0
        sorttri.extend(rows.tolist())

    return sorttri


def be_nice():
    nice(10)
