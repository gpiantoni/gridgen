from nibabel.affines import apply_affine
from numpy import array, ndindex, mean, unravel_index, where, zeros, nansum
from numpy.linalg import norm


def vascular_model(grid, mri, offset):
    n_rows, n_cols = grid.shape

    nd = array(list(ndindex(mri.shape)))
    ndi = apply_affine(mri.affine, nd)

    positions = grid['pos'].reshape(-1, 3)
    chan_xyz = positions + offset

    center_grid = mean(chan_xyz, axis=0)
    dist_chan = norm(ndi - center_grid, axis=1)
    max_dist = max(n_rows, n_cols) * 3 + 5
    i_closeby = dist_chan < max_dist
    good_ndi = ndi[i_closeby]
    val = []
    for pos in chan_xyz:
        dist_chan = norm(good_ndi - pos, axis=1)
        idx = unravel_index(where(i_closeby)[0][dist_chan < 5], mri.shape)
        m = zeros(mri.shape, dtype=int)
        m[idx] = 1
        mq = m * mri.get_data()
        val.append(nansum(mq) / m.sum())

    return array(val).reshape(n_rows, n_cols)
