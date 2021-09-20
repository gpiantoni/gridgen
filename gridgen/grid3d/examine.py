"""Functions to calculate distances and angles between electrodes, to see if
the computed values are within the range of the physical grid."""
from logging import getLogger
from numpy.linalg import norm
from numpy import array, mean, max, min, std, dot, pi, arccos

lg = getLogger(__name__)


def measure_distances(grid2d):
    """Measure the distances between rows and columns (mean, std and range).
    Values are shown in the log

    Parameters
    ----------
    grid2d : (n_rows, n_columns) array
        grid where electrode positions have been computed
    """
    n_rows, n_columns = grid2d.shape

    dist = []
    for i_row in range(n_rows):
        for i_col in range(n_columns - 1):
            dist.append(norm(grid2d[i_row, i_col]['pos'] - grid2d[i_row, i_col + 1]['pos']))
    dist_col = array(dist)

    dist = []
    for i_col in range(n_columns):
        for i_row in range(n_rows - 1):
            dist.append(norm(grid2d[i_row, i_col]['pos'] - grid2d[i_row + 1, i_col]['pos']))
    dist_row = array(dist)

    lg.info(f'Distance between columns: {mean(dist_col):0.3f}mm (sd {std(dist_col):0.3f}mm) [{min(dist_col):0.3f}mm - {max(dist_col):0.3f}mm]')
    lg.info(f'Distance between rows:    {mean(dist_row):0.3f}mm (sd {std(dist_row):0.3f}mm) [{min(dist_row):0.3f}mm - {max(dist_row):0.3f}mm]')


def measure_angles(grid2d):
    """Measure the angles between rows and columns (mean, std and range).
    Values are shown in the log

    Parameters
    ----------
    grid2d : (n_rows, n_columns) array
        grid where electrode positions have been computed. Does not require 'norm', only 'pos'
    """
    n_rows, n_columns = grid2d.shape

    angle = []
    for i_row in range(n_rows - 1):
        for i_col in range(n_columns - 1):
            pos0 = grid2d[i_row, i_col]['pos']
            pos1 = grid2d[i_row + 1, i_col]['pos']
            pos2 = grid2d[i_row, i_col + 1]['pos']
            angle.append(compute_angle(pos0, pos1, pos2))

    for i_row in range(n_rows - 1, 0, -1):
        for i_col in range(n_columns - 1):
            pos0 = grid2d[i_row, i_col]['pos']
            pos1 = grid2d[i_row - 1, i_col]['pos']
            pos2 = grid2d[i_row, i_col + 1]['pos']
            angle.append(compute_angle(pos0, pos1, pos2))

    for i_row in range(n_rows - 1):
        for i_col in range(n_columns - 1, 0, -1):
            pos0 = grid2d[i_row, i_col]['pos']
            pos1 = grid2d[i_row + 1, i_col]['pos']
            pos2 = grid2d[i_row, i_col - 1]['pos']
            angle.append(compute_angle(pos0, pos1, pos2))

    for i_row in range(n_rows - 1, 0, -1):
        for i_col in range(n_columns - 1, 0, -1):
            pos0 = grid2d[i_row, i_col]['pos']
            pos1 = grid2d[i_row - 1, i_col]['pos']
            pos2 = grid2d[i_row, i_col - 1]['pos']
            angle.append(compute_angle(pos0, pos1, pos2))

    lg.info(f'Angle: {mean(angle):0.3f}° (sd {std(angle):0.3f}°) [{min(angle):0.3f}° - {max(angle):0.3f}°]')


_n = lambda x: x / norm(x)


def compute_angle(pos0, pos1, pos2):
    """Compute angle between 3 points. pos0 is the point with the angle of interest

    Parameters
    ----------
    pos0 : (3,) array
        first point (containing the angle)
    pos1 : (3,) array
        second point
    pos2 : (3,) array
        third point

    Returns
    -------
    float
        angle between pos1 and pos2
    """
    return arccos(dot(_n(pos1 - pos0), _n(pos2 - pos0))) / pi * 180
