from logging import getLogger
from numpy.linalg import norm
from numpy import array, mean, max, min, std, dot, pi, arccos

lg = getLogger(__name__)


def measure_distances(grid2d):
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


n = lambda x: x / norm(x)

def compute_angle(pos0, pos1, pos2):
    return arccos(dot(n(pos1 - pos0), n(pos2 - pos0))) / pi * 180
