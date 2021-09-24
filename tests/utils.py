from numpy import array, outer, zeros
from numpy.random import random, seed
from wonambi.utils.simulate import create_data
from nibabel import load, Nifti1Image


def simulate_data():
    chan_name = [f'chan{x + 1}' for x in range(20)]
    data = create_data(chan_name=chan_name, s_freq=1024, signal='sine', sine_freq=70, time=(0, 10))
    ampl_per_chan = array([1, 2, 5, 2, 1, 2, 3, 3, 4, 5, 1, 2, 5, 2, 1, 2, 3, 3, 4, 5])
    data.data[0] = outer(ampl_per_chan, data.data[0][0, :])  # keep same phase
    return data


def simulate_tmap(T1_FILE, TMAP_FILE):

    nii = load(T1_FILE)
    brain = nii.get_fdata()

    roi = zeros(nii.shape, dtype='bool')
    roi[158:, :91, 130:] = 1

    seed(10)
    dat = random(nii.shape).astype('float32') * 100
    dat[(brain <= 0)] = 0
    dat[~roi] = 0

    nifti = Nifti1Image(dat, nii.affine)
    nifti.header.set_xyzt_units('mm')

    nifti.to_filename(str(TMAP_FILE))
