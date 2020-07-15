"""Functions directly useful to compat.py
"""
from numpy import arange, meshgrid, c_, zeros, prod


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
