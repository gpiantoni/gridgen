from numpy import where, array, argmin, zeros, max, c_, ravel
from numpy.linalg import norm
from multiprocessing import Pool
from functools import partial
from scipy.stats import zscore

from nibabel import load
from nibabel.affines import apply_affine


def calculateAngioMap(subj_info, Tthreshold, voxelDepth, plotAngio=False, cortex=None):
    """cortex is not in matlab, but it's necessary"""

    zscore_threshold = 0.1
    xyztCortex, t_surf = voxplot_func_gm(subj_info['sfile'], subj_info['tfile'], cortex, Tthreshold, voxelDepth)
    angioMap = ctmr_vox_plot(cortex, xyztCortex, t_surf, 1, v=None, noplot=True)
    normAngio = (zscore(angioMap) <= zscore_threshold).astype('float')
    return angioMap, normAngio


def voxplot_func_gm(sName, tName, cname, Tthreshold, Dthreshold):

    s_info = load(sName)
    t_info = load(tName)

    s = s_info.get_fdata()
    t = t_info.get_fdata()

    i_t = t >= Tthreshold
    xyz = where_matlab(i_t)
    xyzt = apply_affine(t_info.affine, xyz)

    xyz = where_matlab(s == 1)
    xyzs = apply_affine(s_info.affine, xyz)

    with Pool() as p:
        tsel = p.map(
            partial(
                close_to_surface,
                xyzt=xyzt,
                xyzs=xyzs,
                VoxelDepth=Dthreshold,
            ),
            range(xyzt.shape[0]))
    tsel = array(tsel)

    xyzt = xyzt[tsel]

    with Pool() as p:
        xyztCortex = p.map(
            partial(
                find_closest_vertex,
                cortexpos=cname['pos'],
            ),
            xyzt)

    xyztCortex = array(xyztCortex)

    # matlab order
    t_F = ravel(t, order='F')
    i_t_F = ravel(i_t, order='F')

    t_surf = t_F[i_t_F]
    t_surf = t_surf[tsel]

    return xyztCortex, t_surf


def ctmr_vox_plot(cname, xyz, weights, ssize, v=None, noplot=True):
    """I don't understand implementation but it mirrors the matlab
    implementation.
    """
    cortex = cname

    c = zeros(cortex['pos'].shape[0])
    eps = 1e-5  # we need epsilon for some rounding errors
    for pos, weight in zip(xyz, weights):
        d = (abs(pos - cortex['pos']) <= (ssize + eps)).all(axis=1)
        c = max(c_[c[:, None], d[:, None] * weight], axis=1)

    return c


def close_to_surface(i, xyzt, xyzs, VoxelDepth):
    return (norm(xyzt[i, :] - xyzs, axis=1) <= VoxelDepth).any()


def find_closest_vertex(pos, cortexpos):
    i_min = argmin(norm(pos - cortexpos, axis=1))
    return cortexpos[i_min, :]


def where_matlab(i):
    """WHERE but using matlab convention, in which the last column is sorted first.

    It only works for 3 dimensions
    """
    a = array(where(i)).T
    a = a[a[:, 0].argsort()]
    a = a[a[:, 1].argsort(kind='mergesort')]
    return a[a[:, 2].argsort(kind='mergesort')]
