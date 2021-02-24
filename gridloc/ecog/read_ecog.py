"""Read ECoG data from the grid and convert them to 2D"""
from wonambi import Dataset
from numpy import zeros, NaN, dtype
from wonambi.trans import math, filter_, montage, timefrequency, select
from logging import getLogger
from pathlib import Path

lg = getLogger(__name__)


def read_ecog(ecog_file, begtime=None, endtime=None, bad_channels=None, freq_range=(55, 95)):
    """Read ECoG data

    Parameters
    ----------
    ecog_file : path to file
        Path to ECoG file
    begtime : float
        if specified, data will be read after this point in time (in s)
    endtime : float
        if specified, data will be read until this point in time (in s)
    bad_channels : list of str
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
    lg.debug(f'Reading {ecog_file} between {begtime}s and {endtime}s')
    d = Dataset(Path(ecog_file).resolve())
    data = d.read_data(begtime=begtime, endtime=endtime)

    if bad_channels is not None:
        lg.debug(f'Excluding {len(bad_channels)} channels')
        data = select(data, chan=bad_channels, invert=True)

    lg.debug('Rereference to average')
    data = montage(data, ref_to_avg=True)

    lg.debug('Apply notch filter')
    data = filter_(data, ftype='notch')

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
    ndarray of shape (n_rows, n_columns) with fields:
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
