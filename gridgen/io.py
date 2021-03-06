"""Functions to input and output data
"""
from numpy import (
    abs,
    array,
    cross,
    empty,
    Inf,
    isnan,
    loadtxt,
    NaN,
    savetxt,
    unique,
    where,
    zeros,
    )
from numpy.linalg import norm
from multiprocessing import Pool
from functools import partial
from pathlib import Path
from textwrap import dedent
from warnings import warn

from nibabel import load as nload
from nibabel.freesurfer import read_geometry
from nibabel.affines import apply_affine

from logging import getLogger

from .grid2d import make_grid
from .utils import DTYPE_ECOG, DTYPE_VOLUME

SLICER_HEADER = """# Markups fiducial file version = 4.10
# CoordinateSystem = 0
# columns = id,x,y,z,ow,ox,oy,oz,vis,sel,lock,label,desc,associatedNodeID
"""

lg = getLogger(__name__)


def write_grid2d(grid_file, grid2d):
    """Write the labels of a 2D grid to file

    Parameters
    ----------
    grid_file : Path
        file to write to (preferred extension .csv)
    grid2d : 2d ndarray
        grid (n_rows, n_columns) with fields (label, pos, norm, done)
    """
    with grid_file.open('w') as f:
        for row in grid2d['label']:
            f.write('\t'.join(row) + '\n')


def read_grid2d(grid_file):
    """Read the labels of a 2D grid to file

    Parameters
    ----------
    grid_file : Path
        file to write to (preferred extension .csv)

    Returns
    -------
    grid2d : 2d ndarray
        grid (n_rows, n_columns) with fields (label, pos, norm, done)
    """
    labels = []
    with grid_file.open('r') as f:
        for row in f.readlines():
            labels.append([x.strip() for x in row.split('\t')])

    labels = array(labels)
    grid2d = make_grid(labels.shape[0], labels.shape[1])
    grid2d['label'] = labels
    return grid2d


def write_ecog2d(ecog_file, ecog2d):
    """Write the values of ECoG analysis to file

    Parameters
    ----------
    ecog_file : Path
        file to write to (preferred extension .csv)
    ecog2d : 2d ndarray
        ecog (n_rows, n_columns) with fields (label, ecog, good)
    """
    savetxt(ecog_file, ecog2d['value'], fmt='%.8f', delimiter='\t')


def read_ecog2d(ecog_file, grid_file):
    """Read the values of ECoG analysis

    Parameters
    ----------
    ecog_file : Path
        file with ecog data (in 2d)
    grid_file : Path
        file with labels (in 2d)

    Returns
    -------
    ecog2d : 2d ndarray
        ecog (n_rows, n_columns) with fields (label, ecog)
    """
    ecog = loadtxt(ecog_file, delimiter='\t')

    ecog_on_grid = zeros(ecog.shape, dtype=DTYPE_ECOG)
    ecog_on_grid['value'] = ecog
    ecog_on_grid['good'] = ~isnan(ecog)
    ecog_on_grid['label'] = read_grid2d(grid_file)['label']

    return ecog_on_grid


def read_surf(surf_file, ras_shift=None, normals=True):
    """Read surface file from freesurfer and compute normals.

    Parameters
    ----------
    surf_file : path
        path to Freesurfer pial file (like lh.pial or rh_smooth.pial)
    ras_shift : (3, ) array
        difference in coordinate system between meshes and MRI volume
    normals : bool
        whether to compute normals (it takes one minute roughly)

    Returns
    -------
    dict with fields
        pos : position of each vertex
        tri : triangles (faces) connection between vertices
        pos_norm : normal for each vertex
        tri_norm : normal for each triangle (face)

    Notes
    -----
    All the normals are set to 1, by definition. This is necessary when computing
    the cross-product in later stages.

    It can handle meshes which have vertices that do not belong to a triangle.
    The normal for a vertex without triangle is NaN

    It's strongly adviced to use `ras_shift`, so that all the coordinates are
    consistent in MRI volume space.
    """
    pos, tri = read_geometry(surf_file)
    if ras_shift is None:
        lg.warning('`ras_shift` was not passed. Coordinates in the original mesh coordinate systems, not in the MRI volume space')
    else:
        pos += ras_shift

    surf = {
        'tri': tri,
        'pos': pos,
        'tri_norm': None,
        'pos_norm': None,
        }

    if not normals:
        return surf

    tris = surf['pos'][surf['tri']]

    surf['tri_norm'] = cross(tris[:, 1, :] - tris[:, 0, :], tris[:, 2, :] - tris[:, 0, :])
    surf['tri_norm'] /= norm(surf['tri_norm'], axis=1)[:, None]

    i_vertices = unique(surf['tri'])
    with Pool() as p:
        f_compute = partial(_average_normal_per_vertex, tri_norm=surf['tri_norm'], tri=surf['tri'])
        vert_norm = p.map(f_compute, i_vertices)

    pos_norm = array(vert_norm)
    pos_norm /= norm(pos_norm, axis=1)[:, None]

    surf['pos_norm'] = empty(surf['pos'].shape)
    surf['pos_norm'].fill(NaN)
    surf['pos_norm'][i_vertices] = pos_norm

    return surf


def _average_normal_per_vertex(i, tri_norm, tri):
    return tri_norm[(tri == i).any(axis=1)].mean(axis=0)


def export_grid(grid, ras_shift, grid_file, format=None):
    """
    Parameters
    ----------
    grid : NxNx2x3 array
        grid with positions and normals
    ras_shift : (3, ) array
        shift between MRI volume and meshes
    grid_file : str
        file name to export to (extension is based on format)
    format : str
        'slicer' or 'freeview'. If not specified, both of them

    NOTES
    -----
    Note that the positions should be in the coordinate system of the meshes (not
    the coordinate system of the volume MRI)

    TODO
    ----
    There is something wrong with normals (I guess it depends on how Slicer
    imports them)
    """
    if format is None:
        for one_format in ('freeview', 'slicer'):
            export_grid(grid, ras_shift, grid_file, one_format)
        return

    labels = grid['label'].reshape(-1)
    positions = grid['pos'].reshape(-1, 3) - ras_shift
    normals = grid['norm'].reshape(-1, 3)

    grid_file = Path(grid_file)

    if format == 'slicer':
        grid_file = grid_file.with_suffix('.fcsv')

        with grid_file.open('w') as f:
            f.write(SLICER_HEADER)
            for i in range(positions.shape[0]):
                f.write(f'vtkMRMLMarkupsFiducialNode_{i:03d},{positions[i, 0]:.3f},{positions[i, 1]:.3f},{positions[i, 2]:.3f},{normals[i, 0]:.3f},{normals[i, 1]:.3f},{normals[i, 2]:.3f},1.000,1,0,1,{labels[i]},,\n')

    elif format == 'freeview':
        grid_file = grid_file.with_suffix('.label')

        with grid_file.open('w') as f:
            f.write('#!ascii label  , from subject  vox2ras=TkReg\n')
            f.write(f'{positions.shape[0]:d}\n')
            for i in range(positions.shape[0]):
                f.write(f'{-1:d}  {positions[i, 0]:.3f}  {positions[i, 1]:.3f}  { positions[i, 2]:.3f} 1.000\n')
        warn('make sure that you select the correct volume (.mgz file) associated with the pial surface in freeview when loading the points')


def write_tsv(labels, positions, elec_file):
    """Write electrode position to tsv

    Parameters
    ----------
    labels : list of str
        electrode labels
    positions : NxNx3 or Nx3 array
        electrode position. Make sure to including ras_shift or not if it's in
        MRI space or in mesh space.
    elec_file : Path
        path to write to
    """
    labels = labels.reshape(-1, order='F')
    positions = positions.reshape(-1, 3, order='F')

    elec_file = elec_file.with_suffix('.tsv')
    with elec_file.open('w') as f:
        f.write('name\tx\ty\tz\n')
        for i in range(labels.shape[0]):
            f.write(f'{labels[i]}\t{positions[i, 0]:.3f}\t{positions[i, 1]:.3f}\t{positions[i, 2]:.3f}\n')


def read_surface_ras_shift(T1_path):
    """Freesurfer uses two coordinate systems: one for volumes ("RAS") and
    one for surfaces ("tkReg", "tkRAS", and "Surface RAS").
    To get from surface to volume coordinates, add these numbers.
    To get from volume to surface coordinates, substract these numbers.

    Parameters
    ----------
    T1_path : str
        file path to any .mgz file in the freesurfer folder (e.g. brain.mgz)

    Returns
    -------
    3 array
        offset between volume and surface (tkRAS)
    """
    T1 = nload(str(T1_path))

    return T1.header['Pxyz_c']


def export_transform(offset, transform_file, format='slicer'):
    """Export tkRAS transformation to a transform file.

    Parameters
    ----------
    offset : 3 array
        offset between surface and volume
    transform_file : str
        file name to export to (extension is based on format)
    format : str
        'slicer' or 'freeview'
    """
    assert format == 'slicer'

    transform_file = Path(transform_file)

    transform_file = transform_file.with_suffix('.tfm')

    # use ITK convertion ("ITK's convention is to use LPS coordinate system as opposed to RAS coordinate system in Slicer")
    offset = offset * [-1, -1, 1]

    tfm = """\
        #Insight Transform File V1.0
        #Transform 0
        Transform: AffineTransform_double_3_3
        Parameters: 1 0 0 0 1 0 0 0 1 {:.3f} {:.3f} {:.3f}
        FixedParameters: 0 0 0
        """.format(*offset)

    with transform_file.open('w') as f:
        f.write(dedent(tfm))


def read_volume(volume_file, threshold=-Inf):
    """Read 3D volume. You can also apply threshold to the data

    Parameters
    ----------
    volume_file : path
        path to the nifti file
    threshold : float
        only values above this threshold will be included

    Returns
    -------
    ndarray of shape (n_points, ) with fields:
        - pos : 3 floats (specifying the x, y, z position)
        - value : float (actual value at that point)

    Notes
    -----
    Values below the threhold are not included in the output
    """
    volume = nload(str(volume_file))
    data = volume.get_fdata()
    if threshold is not None:
        i = data >= threshold
    else:
        i = abs(data) > 0.001  # exclude small values

    output = zeros(i.sum(), dtype=DTYPE_VOLUME)

    output['pos'] = apply_affine(volume.affine, array(where(i)).T)
    if threshold is not None:
        output['value'] = 1
    else:
        output['value'] = data[i]

    return output


def read_mri(T1_file, dura_file, pial_file=None, func_file=None, func_threshold=None):
    """
    """
    ras_shift = read_surface_ras_shift(T1_file)
    lg.debug(f'Reading positions and computing normals of {dura_file}')
    dura = read_surf(dura_file, ras_shift=ras_shift)

    if pial_file is None:
        pial = None
    else:
        lg.debug(f'Reading positions of {pial_file}')
        pial = read_surf(pial_file, normals=False, ras_shift=ras_shift)

    if func_file is None:
        func = None
    else:
        lg.debug(f'Reading functional image from {func_file} and thresholding at {func_threshold}')
        func = read_volume(func_file, func_threshold)
        if func_threshold is not None:
            func['value'] = 1  # set the value of each voxel above threshold to 1

    out = {
        "ras_shift": ras_shift,
        "dura": dura,
        "pial": pial,
        "func": func,
        }
    return out


def export_electrodes(output, model, mris):
    grid_file = output / 'electrodes'
    export_grid(model['grid'], mris['ras_shift'], grid_file)
    write_tsv(model['grid']['label'], model['grid']['pos'], grid_file)
    lg.debug(f'Exported electrodes to {grid_file} (coordinates in MRI volume space, not mesh space)')
