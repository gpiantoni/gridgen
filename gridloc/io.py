from numpy import cross, array
from numpy.linalg import norm
from multiprocessing import Pool
from functools import partial
from pathlib import Path

from nibabel.freesurfer import read_geometry

SLICER_HEADER = """# Markups fiducial file version = 4.10
# CoordinateSystem = 0
# columns = id,x,y,z,ow,ox,oy,oz,vis,sel,lock,label,desc,associatedNodeID
"""


def read_surf(surf_file, normals=True, norm_to_one=True):

    pos, tri = read_geometry(surf_file)
    surf = {
        'tri': tri,
        'pos': pos,
        'tri_norm': None,
        'pos_norm': None,
        }

    tris = surf['pos'][surf['tri']]

    if not normals:
        return surf

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


def export_grid(grid, grid_file, format='slicer'):
    """
    Parameters
    ----------
    grid : NxNx2x3 array
        grid with positions and normals
    grid_file : str
        file name to export to (extension is based on format)
    format : str
        'slicer' or 'freeview'

    TODO
    ----
    There is something wrong with normals (I guess it depends on how Slicer
    imports them)
    """
    positions = grid[:, :, 0, :].reshape(-1, 3)
    normals = grid[:, :, 1, :].reshape(-1, 3)

    grid_file = Path(grid_file)

    if format == 'slicer':
        grid_file = grid_file.with_suffix('.fcsv')

        with grid_file.open('w') as f:
            f.write(SLICER_HEADER)
            for i in range(positions.shape[0]):
                f.write(f'vtkMRMLMarkupsFiducialNode_{i:03d},{positions[i, 0]:.3f},{positions[i, 1]:.3f},{positions[i, 2]:.3f},{normals[i, 0]:.3f},{normals[i, 1]:.3f},{normals[i, 2]:.3f},1.000,1,1,1,ELEC{i + 1:03d},,\n')

    elif format == 'freeview':
        grid_file = grid_file.with_suffix('.label')

        with grid_file.open('w') as f:
            f.write('#!ascii label  , from subject  vox2ras=TkReg\n')
            f.write(f'{positions.shape[0]:d}\n')
            for i in range(positions.shape[0]):
                f.write(f'{-1:d}  {positions[i, 0]:.3f}  {positions[i, 1]:.3f}  { positions[i, 2]:.3f} 1.000\n')
        print('make sure that you select the brain.mgz associated with the pial surface in freeview')
