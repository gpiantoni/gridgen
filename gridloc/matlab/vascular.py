from numpy import where, array, argmin
from numpy.linalg import norm
from multiprocessing import Pool
from functools import partial

from nibabel import load
from nibabel.affines import apply_affine


def voxplot_func_gm(sName, tName, cname, Tthreshold, Dthreshold):

    s_info = load(sName)
    t_info = load(tName)

    s = s_info.get_fdata()
    t = t_info.get_fdata()

    i_t = t >= Tthreshold
    t_surf = t[i_t]
    xyz = where(i_t)
    xyz = array(xyz).T
    xyzt = apply_affine(t_info.affine, xyz)

    xyz = where(s == 1)
    xyz = array(xyz).T
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

    tsel = array([x for x in tsel if x is not None])
    xyzt = xyzt[tsel]
    t_surf = t_surf[tsel]

    with Pool() as p:
        xyztCortex = p.map(
            partial(
                find_closest_vertex,
                cortexpos=cname['pos'],
            ),
            xyzt)

    xyztCortex = array(xyztCortex)

    return xyztCortex, t_surf


def close_to_surface(i, xyzt, xyzs, VoxelDepth):
    if (norm(xyzt[i, :] - xyzs, axis=1) <= VoxelDepth).any():
        return i

def find_closest_vertex(pos, cortexpos):
    i_min = argmin(norm(pos - cortexpos, axis=1))
    return cortexpos[i_min, :]
