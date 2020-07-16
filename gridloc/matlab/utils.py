"""Functions directly useful to compat.py
"""

from numpy import arange, meshgrid, c_, zeros, prod, argmax, dot, cross, NaN, array
from numpy.linalg import norm
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
    """Taken from plane_intersect.m"""
    P = array([0., 0., 0.])
    N = cross(N1, N2)
    if norm(N) < EPSILON:  # parallel or coincide
        return NaN, NaN

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
