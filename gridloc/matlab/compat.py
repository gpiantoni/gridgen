"""The functions in this module should have the same name and argument signature
of the Matlab functions, for compatibility"""

from numpy import array
from numpy.linalg import norm
from .geometry import project_to_cortex


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

    Notes
    -----
    It's not possible to replicate projectElectrodes.m completely. When `normway`
    is a scalar, it takes the normal of the points of the mesh which are closer
    than `normway`. However, some points have a normal of (0, 0, 0) (default assigned
    if the vertex does not belong to any triangle). projectElectrodes.m includes
    those (0, 0, 0) in the calculation, but it might not be correct.
    See l. 138 (there are no NaN in normals but only (0, 0, 0).

    To replicate it completely, you'd need to take the empty normals into account.
    """
    normdist = normway  # if normway is scalar

    normals = []
    pints = []
    for electrode in subjstructs['electrodes']:
        dvect = norm(electrode - surf['pos'], axis=1)  # l. 104-112
        closevert = dvect < (normdist ** 2)  # l. 120
        normal = surf['pos_norm'][closevert, :].mean(axis=0)  # l. 144  (but see Notes)
        pint = project_to_cortex(surf, electrode, normal)[1]  # l. 198-237

        normals.append(normal)
        pints.append(pint)

    subjstructs['normals'] = array(normals)
    subjstructs['trielectrodes'] = array(pints)

    return subjstructs
