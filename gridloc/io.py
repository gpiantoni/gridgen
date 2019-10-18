from numpy import cross, array
from numpy.linalg import norm
from multiprocessing import Pool
from functools import partial

from nibabel.freesurfer import read_geometry

def read_surf(surf_file, norm_to_one=False):

    pos, tri = read_geometry(surf_file)
    surf = {
        'tri': tri,
        'pos': pos,
        }

    tris = surf['pos'][surf['tri']]
    surf['tri_norm'] = cross(tris[:, 1, :] - tris[:, 0, :], tris[:, 2, :] - tris[:, 0, :])

    if norm_to_one:
        surf['tri_norm'] /= norm(surf['tri_norm'], axis=1)[:, None]

    with Pool() as p:
        f_compute = partial(_average_normal_per_vertex, tri_norm=surf['tri_norm'], tri=surf['tri'])
        vert_norm = p.map(f_compute, range(surf['pos'].shape[0]))

    surf['pos_norm'] = array(vert_norm)

    if norm_to_one:
        surf['pos_norm'] /= norm(surf['pos_norm'], axis=1)[:, None]

    return surf


def _average_normal_per_vertex(i, tri_norm, tri):
    return tri_norm[(tri == i).any(axis=1)].mean(axis=0)
