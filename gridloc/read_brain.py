from wonambi import Dataset
from wonambi.trans.select import _create_subepochs
from numpy import std, median, zeros, array
from wonambi.trans import math, filter_, montage, timefrequency, select


def read_brain(ecog_file, begtime=None, endtime=None, freq_range=(55, 95),
               chan_std_threshold=100, grid_dims=(16, 8)):

    d = Dataset(ecog_file)
    data = d.read_data(begtime=begtime, endtime=endtime)
    dd = montage(data, ref_chan=['chan1', 'chan128'])
    filt = filter_(dd, ftype='notch')
    x = _create_subepochs(
        data.data[0],
        int(data.s_freq),
        int(data.s_freq / 2))
    s = median(std(x, axis=2), axis=1)

    good_chans = array(d.header['chan_name'])[s < chan_std_threshold]
    print(len(good_chans))

    filtx = select(filt, chan=good_chans)

    f = timefrequency(filtx, method='spectrogram', duration=2, overlap=0.5, taper='hann')
    ff = math(f, operator_name='median', axis='time')
    ff = math(ff, operator_name='dB')
    fff = select(ff, freq=freq_range)
    x = math(fff, operator_name='mean', axis='freq')

    i = 0
    gamma = zeros(grid_dims)
    for i_y in range(grid_dims[1]):
        for i_x in range(grid_dims[0]):
            i += 1
            gamma[i_x, i_y] = x(chan=f'chan{i}', trial=0)

    return gamma
