from wonambi import Dataset
from wonambi.trans.select import _create_subepochs
from numpy import std, median, zeros, array, NaN
from wonambi.trans import math, filter_, montage, timefrequency, select
from logging import getLogger

lg = getLogger(__name__)


def read_brain(ecog_file, begtime=None, endtime=None, ref_chan=None,
               freq_range=(55, 95), chan_std_threshold=100, grid_dims=(16, 8)):

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
    lg.info(f'Keeping {len(good_chans)} channels')

    data = select(data, chan=good_chans)

    tf = timefrequency(data, method='spectrogram', duration=2, overlap=0.5, taper='hann')
    tf = math(tf, operator_name='median', axis='time')
    tf = math(tf, operator_name='dB')
    tf = select(tf, freq=freq_range)
    tf = math(tf, operator_name='mean', axis='freq')

    return tf

"""
    i = 0
    gamma = zeros(grid_dims)
    for i_y in range(grid_dims[1]):
        for i_x in range(grid_dims[0]):
            i += 1
            gamma[i_x, i_y] = x(chan=f'chan{i}', trial=0)

    # TODO: clean up this part
    gamma[0, 0] = NaN
    gamma[-1, -1] = NaN
"""


def main(parameters):

    gamma = read_brain(
        parameters['ecog']['file'],
        parameters['ecog']['begtime'],
        parameters['ecog']['endtime'],
        parameters['ecog']['ref_chan'],
        parameters['ecog']['freq_range'],
        parameters['ecog']['chan_std_threshold'],
        (parameters['grid']['n_rows'], parameters['grid']['n_columns']))

    return gamma
