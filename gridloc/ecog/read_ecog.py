from wonambi import Dataset
from wonambi.trans.select import _create_subepochs
from numpy import std, median, zeros, array, NaN, dtype, setdiff1d
from wonambi.trans import math, filter_, montage, timefrequency, select
from logging import getLogger

lg = getLogger(__name__)


def read_brain(ecog_file, begtime=None, endtime=None, ref_chan=None,
               freq_range=(55, 95), chan_std_threshold=100):

    lg.debug(f'Reading {ecog_file} between {begtime}s and {endtime}s')
    d = Dataset(ecog_file)
    data = d.read_data(begtime=begtime, endtime=endtime)

    if ref_chan is not None:
        lg.debug(f'Rereference to {ref_chan}')
        data = montage(data, ref_chan=ref_chan)

    lg.debug('Apply notch filter')
    data = filter_(data, ftype='notch')

    lg.debug(f'Computing the median of the s.d. in subepochs')
    x = _create_subepochs(
        data.data[0],
        int(data.s_freq),
        int(data.s_freq / 2))
    s = median(std(x, axis=2), axis=1)

    good_chans = array(d.header['chan_name'])[s < chan_std_threshold]
    lg.debug(f'Rejecting {len(d.header["chan_name"]) - len(good_chans)} channels')
    good_chans = setdiff1d(good_chans, ref_chan)
    lg.debug(f'Removing {len(ref_chan)} channels')
    lg.info(f'Using {len(good_chans)} channels')

    data = select(data, chan=good_chans)

    lg.debug(f'Computing Spectrogram')
    tf = timefrequency(data, method='spectrogram', duration=2, overlap=0.5, taper='hann')
    tf = math(tf, operator_name='median', axis='time')
    tf = math(tf, operator_name='dB')
    tf = select(tf, freq=freq_range)
    tf = math(tf, operator_name='mean', axis='freq')

    return tf


def put_ecog_on_grid2d(ecog, grid2d):
    d_ = dtype([
        ('label', '<U256'),
        ('ecog', 'f4'),
        ('good', 'bool'),
        ])
    ecog_on_grid = zeros(grid2d.shape, dtype=d_)

    ecog_on_grid['label'] = grid2d['label']
    ecog_on_grid['ecog'].fill(NaN)

    for x in range(ecog_on_grid.shape[0]):
        for y in range(ecog_on_grid.shape[1]):
            label = ecog_on_grid['label'][x, y]
            if label in ecog.axis['chan'][0]:
                ecog_on_grid['ecog'][x, y] = ecog(trial=0, chan=label)
                ecog_on_grid['good'][x, y] = True

    return ecog_on_grid
