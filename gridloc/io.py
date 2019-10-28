from numpy import cross, array
from numpy.linalg import norm
from multiprocessing import Pool
from functools import partial

from nibabel.freesurfer import read_geometry

SLICER_HEADER = """# Markups fiducial file version = 4.10
# CoordinateSystem = 0
# columns = id,x,y,z,ow,ox,oy,oz,vis,sel,lock,label,desc,associatedNodeID
"""


def read_surf(surf_file, norm_to_one=True):

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


def export_grid_to_3dslicer(grid, grid_file):
    """
    TODO
    ----
    There is something wrong with normals (I guess it depends on how Slicer
    imports them
    """
    positions = grid[:, :, 0, :].reshape(-1, 3)
    normals = grid[:, :, 1, :].reshape(-1, 3)

    with open(grid_file, 'w') as f:
        f.write(SLICER_HEADER)
        for i in range(positions.shape[0]):
            f.write(f'vtkMRMLMarkupsFiducialNode_{i:03d},{positions[i, 0]:.3f},{positions[i, 1]:.3f},{positions[i, 2]:.3f},{normals[i, 0]:.3f},{normals[i, 1]:.3f},{normals[i, 2]:.3f},1.000,1,1,1,ELEC{i + 1:03d},,\n')
