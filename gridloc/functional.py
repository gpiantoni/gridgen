from functools import partial
from itertools import product
from logging import getLogger
from multiprocessing import Pool
import warnings

from numpy import (ndindex,
                   array,
                   sum,
                   power,
                   zeros,
                   NaN,
                   isfinite,
                   nansum,
                   isnan,
                   where,
                   )
from numpy.linalg import norm, inv
from scipy.stats import norm as normdistr
from nibabel.affines import apply_affine
from nibabel import load as nload

lg = getLogger(__name__)

# 1 sigma = 0.6065306597126334


def calc_fmri_at_elec(measure_nii, electrodes_file, distance, kernels, graymatter, output_dir):
    """
    Calculate the (weighted) average of fMRI values at electrode locations
    """

    img = nload(str(measure_nii))
    mri = img.get_data()
    mri[mri == 0] = NaN

    labels = electrodes.electrodes.get(map_lambda=lambda x: x['name'])
    chan_xyz = array(electrodes.get_xyz())

    nd = array(list(ndindex(mri.shape)))
    ndi = from_mrifile_to_chan(img, nd)

    if graymatter:
        gm_mri = nload(str(graymatter)).get_data().astype(bool)
        mri[~gm_mri] = NaN

    lg.debug(f'Computing fMRI values for {measure_nii.name} at {len(labels)} electrodes and {len(kernels)} "{distance}" kernels')
    fmri_vals, n_voxels = compute_kernels(kernels, chan_xyz, mri, ndi, distance)

    fmri_vals_tsv = output_dir / replace_underscore(measure_nii.name, 'compare.tsv')
    n_voxels_tsv = output_dir / replace_underscore(measure_nii.name, 'nvoxels.tsv')

    with fmri_vals_tsv.open('w') as f:
        f.write('channel\t' + '\t'.join(str(one_k) for one_k in kernels) + '\n')
        for one_label, val_at_elec in zip(labels, fmri_vals):
            f.write(one_label + '\t' + '\t'.join(str(one_val) for one_val in val_at_elec) + '\n')

    with n_voxels_tsv.open('w') as f:
        f.write('channel\t' + '\t'.join(str(one_k) for one_k in kernels) + '\n')
        for one_label, val_at_elec in zip(labels, n_voxels):
            f.write(one_label + '\t' + '\t'.join(str(one_val) for one_val in val_at_elec) + '\n')

    return fmri_vals_tsv, n_voxels_tsv


def from_chan_to_mrifile(img, xyz):
    return apply_affine(inv(img.affine), xyz).astype(int)


def from_mrifile_to_chan(img, xyz):
    return apply_affine(img.affine, xyz)


def compute_kernels(kernels, chan_xyz, mri, ndi, distance, n_cpu=None):
    partial_compute_chan = partial(compute_chan, ndi=ndi, mri=mri, distance=distance)

    args = product(chan_xyz, kernels)
    if n_cpu is None:
        fmri_vals = [partial_compute_chan(*arg) for arg in args]
    else:
        lg.debug(f'Number of CPU: {n_cpu}')
        with Pool(n_cpu) as p:
            fmri_vals = p.starmap(partial_compute_chan, args)

    fmri_vals = array(fmri_vals).reshape(-1, len(kernels))

    return fmri_vals


def compute_chan(pos, KERNEL, ndi, mri, distance):
    dist_chan = norm(ndi - pos, axis=1)

    if distance == 'gaussian':
        m = normdistr.pdf(dist_chan, scale=KERNEL)
        m /= normdistr.pdf(0, scale=KERNEL)  # normalize so that peak is at 1, so that it's easier to count voxels

    elif distance == 'sphere':
        m = zeros(dist_chan.shape)
        m[dist_chan <= KERNEL] = 1

    elif distance == 'inverse':
        m = power(dist_chan, -1 * KERNEL)

    m = m.reshape(mri.shape)
    m[isnan(mri)] = NaN

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        m /= sum(m[isfinite(m)])  # normalize so that the sum of the finite numbers is 1

    mq = m * mri
    return nansum(mq)
