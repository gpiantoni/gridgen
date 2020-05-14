from numpy import array, outer
from wonambi.utils.simulate import create_data


def simulate_data():
    chan_name = [f'chan{x + 1}' for x in range(20)]
    data = create_data(chan_name=chan_name, s_freq=1024, signal='sine', sine_freq=70, time=(0, 10))
    ampl_per_chan = array([1, 2, 5, 2, 1, 2, 3, 3, 4, 5, 1, 2, 5, 2, 1, 2, 3, 3, 4, 5])
    data.data[0] = outer(ampl_per_chan, data.data[0][0, :])  # keep same phase
    data.export(ECOG_FILE, export_format='BrainVision')
