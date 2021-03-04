"""The functions in this module should have the same name and argument signature
of the Matlab functions, for compatibility"""

from numpy import array, arange, meshgrid, nanmean, isnan, dot, prod, arctan2, cross, pi, mean, errstate
from numpy.linalg import norm
from multiprocessing import Pool
from functools import partial

from .geometry import project_to_cortex
from .utils import plane_intersect, AxelRot, _apply_affine, _sort_closest_triangles, calcCoords
from ..utils import be_nice

from ..geometry import calc_plane_to_axis
from ..construct import make_grid_with_labels


def getROI(surf, ref_vert, ROIsize=18, intElec=3):
    """generate 10x10 grid, on one plane

    TODO
    ----
    we should calculate pos and norm from
    pos -> projected position
    normal -> normal to the hullcortex
    """
    pos = surf['pos'][ref_vert, :]
    normal = surf['pos_norm'][ref_vert, :]

    # trying to replicate getROI.m (l. 69-90) but we do not consider the orientation here
    steps = arange(-ROIsize + intElec, ROIsize - intElec, intElec)
    x_mesh, y_mesh = meshgrid(steps, steps)

    ROI = []
    for x, y in zip(x_mesh.flatten(), y_mesh.flatten()):
        coords_2d = array([x, y])
        plane = calc_plane_to_axis(normal)
        target = coords_2d @ plane + pos
        ROI.append(target)

    return array(ROI)


def projectElectrodes(surf, subjstructs, normway, normUse=False, interstype='',
                      intersval=0):
    """Replicate projectElectrodes. This function is a ray-triangle intersection.
    The normal of the point is specified in `normway` if normway has three values.
    Otherwise, it computes the normal based on the direction of the points of the
    mesh which are closer than `normway`.

    Parameters
    ----------
    surf : dict with 'pos', 'tri'
        surface of the brain use to project the electrodes (it's not necessary
        to have 'tri_norm')
    subjstructs : dict with 'electrodes'
        location of the electrodes
    normway : float
        distance to use to include mesh points when computing the normal
    normUse : bool
        if true, it uses the normals of `subjstructs`, if false it recomputes them
    interstype : str
        if `''` all triangles of the model are processed; if `'fixed'`, you need
        to specify a radius `intersval`)
    intersval : float
        radius when `interstype` == `'fixed'`

    Returns
    -------
    dict with 'electrodes', 'normals', 'trielectrodes'
        for each electrode, it returns its normal (based on neighboring mesh
        vertices) and the projection onto the surface.
    """
    assert interstype in ('', 'fixed')

    normdist = normway  # if normway is scalar

    normals = []
    pints = []
    for i in range(subjstructs['electrodes'].shape[0]):
        electrode = subjstructs['electrodes'][i, :]

        if normUse:
            normal = subjstructs['normal'][i, :]
        else:
            normal = normElec(surf, electrode, normdist)

        if interstype == 'fixed':
            sorttri = _sort_closest_triangles(surf, electrode, intersval)

        else:
            sorttri = None

        # note the sign of normal is inverted
        pint = project_to_cortex(surf, electrode, normal, sorted_triangles=sorttri)[1]  # l. 198-237

        normals.append(normal)
        pints.append(pint)

    subjstructs['normal'] = array(normals)
    subjstructs['trielectrodes'] = array(pints)

    return subjstructs


def normElec(surf, electrode, normdist, NaN_as_zeros=True):
    """
    Notes
    -----
    When `normway` is a scalar, it takes the normal of the points of the mesh which are closer
    than `normway`. However, some points have a normal of (0, 0, 0) (default assigned
    if the vertex does not belong to any triangle). projectElectrodes.m includes
    those (0, 0, 0) in the calculation, but it might not be correct.
    See l. 138 (there are no NaN in normals but only (0, 0, 0)).

    To replicate the matlab behavior, make sure that `NaN_as_zeros` is True.
    """
    dvect = norm(electrode - surf['pos'], axis=1)  # l. 104-112 of projectElectrodes.m
    closevert = dvect < normdist  # l. 120 of projectElectrodes.m
    normal = surf['pos_norm'][closevert, :].mean(axis=0)  # l. 144 of projectElectrodes.m

    normals2av = surf['pos_norm'][closevert, :].copy()
    if NaN_as_zeros:
        normals2av[isnan(normals2av)] = 0
    normal = nanmean(normals2av, axis=0)

    return normal


def calcTangent(hullcortex, c, coords, dims, lngth, hemi):

    N = normElec(hullcortex, c, 25)
    d = dot(N, c)

    point = c
    normal = N
    tangPlane = coords.copy().reshape(prod(dims), 3)
    tangPlane[:, 0] = (-normal[1] * tangPlane[:, 1] - normal[2] * tangPlane[:, 2] - d) / normal[0]  # l. 19
    N = tangPlane[0, :] - tangPlane[7, :]
    M = coords[0, :] - coords[7, :]

    AngleBetweenPlanes = arctan2(norm(cross(M, N)), dot(M, N))

    Point, IntersectionPlanes = plane_intersect(array([1., 0., 0.]), coords[0, :], normal, point)

    RotMatrix = AxelRot(AngleBetweenPlanes, IntersectionPlanes, Point)
    coords_Tangent = _apply_affine(coords, RotMatrix)

    if hemi == 'r':
        coords_Tangent[:, 0] = coords_Tangent[::-1, 0]  # l. 54-65

    return coords_Tangent


def createGrid(sub, rotation=None, turns=None, auxDims=(8, 16), subj_info=None, hullcortex=None):
    """Create GRID per ROI point and project on cortex

    Parameters
    ----------
    sub : dict of array
        the same as 'projectedROIpoints'
    rotation : None
        not used (but it might be in the future)
    turns : None
        not used (but it might be in the future)
    auxDims : tuple of int
        dimensions of the grid
    subj_info : dict
        information about this subject
    hullcortex : dict with 'pos' and 'tri'
        hull cortex
    """
    hemi = subj_info['hemiVect']['hemi']
    intElec = subj_info['intElec']
    f = partial(createGrid_per_point, sub=sub, intElec=intElec, auxDims=auxDims, hemi=hemi, hullcortex=hullcortex)
    with Pool(initializer=be_nice) as p:
        ROI = p.map(f, range(sub['electrodes'].shape[0]))

    return ROI


def createGrid_per_point(roi_punt, sub, intElec, auxDims, hemi, hullcortex):
    ROI_pos = sub['trielectrodes'][roi_punt, :]
    ROI_norm = sub['normal'][roi_punt, :]

    twoDgridElectrodes = calcCoords(ROI_pos, intElec, auxDims)
    electrodes_start = calcTangent(hullcortex, ROI_pos, twoDgridElectrodes, auxDims, prod(auxDims), hemi)

    normNormals = ROI_norm / norm(ROI_norm)

    # apply first rotation
    _transform = AxelRot(-45 / 180 * pi, normNormals, ROI_pos)
    electrodes_start = _apply_affine(electrodes_start, _transform)  # l. 44-50

    normdist = 25
    intersval = 30  # use largest range anyway

    coords = []
    rotations = range(90)
    for rotation in rotations:
        _transform = AxelRot(rotation / 180 * pi, normNormals, ROI_pos)
        new_turn = {'electrodes': _apply_affine(electrodes_start, _transform)}
        new_turn = projectElectrodes(hullcortex, new_turn, normdist, normUse=False, interstype='fixed', intersval=intersval)
        coords.append(new_turn)

    return coords


def projectToCoarser(ROI, cortexcoarser, turns=None):
    """
    Project the create grid onto the coarser model

    turns is not necessary
    """
    f = partial(projectToCoarser_per_point, cortexcoarser=cortexcoarser)
    with Pool() as p:
        ROI = p.map(f, ROI)

    return ROI


def projectToCoarser_per_point(coords, cortexcoarser):
    intersval = 35
    for one_coord in coords:
        ss = {
            'electrodes': one_coord['trielectrodes'],
            'normal': one_coord['normal'],
            }
        ss = projectElectrodes(cortexcoarser, ss, 25, normUse=True, interstype='fixed', intersval=intersval)
        one_coord['pint'] = ss['trielectrodes']

    return coords


def calculateModel(null, ROI, cortex, normAngio=None):
    """Use only the first dimension
    Cortex is only needed when using normAngio
    """
    for coords in ROI:
        for one_coord in coords:

            # only the first dimension
            McM = abs(one_coord['pint'][:, 0] - one_coord['trielectrodes'][:, 0])
            McM = 1 - (McM - McM.min()) / McM.max()
            one_coord['McM'] = McM

            if normAngio is None:
                one_coord['weights'] = one_coord['McM']

    if normAngio is not None:
        f = partial(_calculateVascularModel, cortex=cortex, normAngio=normAngio)
        with Pool() as p:
            ROI = p.map(f, ROI)

    return ROI


def _calculateVascularModel(coords, cortex, normAngio):
    for one_coord in coords:
        v = []
        for pint in one_coord['pint']:
            dist = norm(pint - cortex['pos'], axis=1)
            with errstate(invalid='ignore'):
                idx = dist <= 2
            v.append(mean(normAngio[idx]))
        one_coord['MvM'] = array(v)

        one_coord['weights'] = 0.5 * one_coord['MvM'] + 0.5 * one_coord['McM']

    return coords


def indexFuncLegacy(subj_info):
    """Computes the label to assign to the electrode position (so not to the
    ecog activity, as in matlab function)

    Requires
      - hemiVect:
        - hemi: 'l' or 'r'
        - side: 'u' (up), 'd' (down), 'l' (left), or 'r' (right)

    - dims [n_columns, n_rows]
    Note that `dims` is first columns, then rows, so for the normal grid it's
    [8, 16]

    """
    n_columns, n_rows = subj_info['dims']

    if subj_info['hemiVect']['hemi'] == 'l':
        if subj_info['hemiVect']['side'] == 'u':
            direction = 'TBLR'

        elif subj_info['hemiVect']['side'] == 'd':
            direction = 'BTRL'

        elif subj_info['hemiVect']['side'] == 'l':
            direction = 'LRBT'
            n_rows, n_columns = n_columns, n_rows

        elif subj_info['hemiVect']['side'] == 'r':
            direction = 'RLTB'
            n_rows, n_columns = n_columns, n_rows

    elif subj_info['hemiVect']['hemi'] == 'r':
        if subj_info['hemiVect']['side'] == 'u':
            direction = 'TBRL'

        elif subj_info['hemiVect']['side'] == 'd':
            direction = 'BTLR'

        elif subj_info['hemiVect']['side'] == 'r':
            direction = 'LRTB'
            n_rows, n_columns = n_columns, n_rows

        elif subj_info['hemiVect']['side'] == 'l':
            direction = 'RLBT'
            n_rows, n_columns = n_columns, n_rows

    elec2d = make_grid_with_labels(n_rows, n_columns, direction=direction, chan_pattern='chan{}')
    return elec2d['label'].flatten(order='F')
