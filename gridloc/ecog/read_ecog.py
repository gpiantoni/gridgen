from wonambi import Dataset
from numpy import zeros, NaN, dtype
from wonambi.trans import math, filter_, montage, timefrequency, select
from logging import getLogger
from pathlib import Path

lg = getLogger(__name__)


def read_ecog(file, begtime=None, endtime=None, ref_chan=None, bad_chan=None,
              freq_range=(55, 95)):
    """Read ECoG data

    Parameters
    ----------
    file : path to file
        Path to ECoG file
    begtime : float
        if specified, data will be read after this point in time (in s)
    endtime : float
        if specified, data will be read until this point in time (in s)
    ref_chan : list of str
        list of channels to use as reference
    bad_chan : list of str
        list of channels to exclude
    freq_range : two floats
        frequency range to compute the average power spectrum

    Returns
    -------
    wonambi Data
        where you get one value per channels, which is the average in the
        frequency range

    TODO
    ----
    check whether we should use 'mean' or 'median' across time
    """
    lg.debug(f'Reading {file} between {begtime}s and {endtime}s')
    d = Dataset(Path(file).resolve())
    data = d.read_data(begtime=begtime, endtime=endtime)

    if ref_chan is not None:
        lg.debug(f'Rereference to {ref_chan}')
        data = montage(data, ref_chan=ref_chan)

    lg.debug('Apply notch filter')
    data = filter_(data, ftype='notch')

    if bad_chan is not None:
        data = select(data, chan=bad_chan, invert=True)

    lg.debug('Computing Spectrogram')
    tf = timefrequency(data, method='spectrogram', duration=2, overlap=0.5, taper='hann')
    tf = math(tf, operator_name='mean', axis='time')
    tf = math(tf, operator_name='dB')

    lg.debug(f'Selecting frequency range {freq_range[0]:02.2f}-{freq_range[1]:02.2f}')
    tf = select(tf, freq=freq_range)
    tf = math(tf, operator_name='mean', axis='freq')

    return tf


def put_ecog_on_grid2d(ecog, grid2d):
    """Arrange values of the ecog in a 2d grid, based on the 2d electrode
    location

    Parameters
    ----------
    ecog : wonambi Data
        output of read_ecog(). There should be one value per channel
    grid2d : 2d array
        2d array with the channel labels (input is a 2d array with one field called 'label')

    Returns
    -------
    2d array
        array with same shape as grid2d, with fields:
        - label : channel label
        - ecog : value computed from 'ecog'
        - good : whether to include the channel or not
    """
    d_ = dtype([
        ('label', '<U256'),
        ('ecog', 'f8'),
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
