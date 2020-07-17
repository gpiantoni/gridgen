"""The functions in this module should have the same name and argument signature
of the Matlab functions, for compatibility"""

from numpy import array, arange, meshgrid, nanmean, isnan, dot, prod, arctan2, cross
from numpy.linalg import norm
from ..geometry import calc_plane_to_axis
from .geometry import project_to_cortex
from .utils import plane_intersect, AxelRot, _apply_affine


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


def projectElectrodes(surf, subjstructs, normway):
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

    Returns
    -------
    dict with 'electrodes', 'normals', 'trielectrodes'
        for each electrode, it returns its normal (based on neighboring mesh
        vertices) and the projection onto the surface.
    """
    normdist = normway  # if normway is scalar

    normals = []
    pints = []
    for electrode in subjstructs['electrodes']:
        normal = normElec(surf, electrode, normdist)
        # note the sign of normal is inverted
        pint = project_to_cortex(surf, electrode, -1 * normal)[1]  # l. 198-237

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
