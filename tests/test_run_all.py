from numpy import argmin
from gridloc.search import calc_plane_to_axis

from gridloc.construct import construct_grid
from gridloc.compare import compare_models

from itertools import product

from multiprocessing import Pool
import mkl
mkl.set_num_threads(2)


def _compute_grid(start_vert, rotation):
    grid = construct_grid(surf, start_vert, start_label, n_rows, n_cols, rotation=rotation)
    return compare_models(grid, pial, angio, offset, gamma)

def search_grid():
    pos = surf['pos'][start_vert, :]
    normal = surf['pos_norm'][start_vert, :]

    vec = []
    delta = 1
    for x in range(4):
        for y in range(4):
            coords_2d = array([x, y]) * delta
            plane = calc_plane_to_axis(normal)
            target = coords_2d @ plane + pos
            new_vector = argmin(norm(surf['pos'] - target, axis=1))
            vec.append(new_vector)


def notest():
    from scipy.io import loadmat
    gamma = gamma.reshape(8, 16).T

    surf = read_surf(surf_file)

    offset = array([-5.5789, 4.33667, 33.9133]) * -1

    pial = read_surf(pial_file, normals=False)

    args = product(vec, list(range(90 - 15, 90 + 15, 3)))
    with Pool() as p:
        c = p.starmap(_compute_grid, args)
    x = array(c)
    x = x.reshape(10, 4, 4, order='F')


n_rows = 16
n_cols = 8
    start_label = 'elec079'

    order = 'minor'
    rotation = 60
    start_vert = 14662 # 0.340
    # start_vert = 38242
    start_vert = 13873
