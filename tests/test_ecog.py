from numpy import array, outer
from wonambi.utils.simulate import create_data

from gridloc.ecog.read_ecog import read_ecog, put_ecog_on_grid2d
from gridloc.construct import make_grid_with_labels

from .paths import ECOG_FILE


def test_read_ecog():
    simulate_data()
    tf = read_ecog(ECOG_FILE)
    assert tf.data[0][0] < tf.data[0][1] < tf.data[0][2]

    grid2d = make_grid_with_labels(4, 5, 'TBLR', 'chan{}')
    ecog_grid = put_ecog_on_grid2d(tf, grid2d)

    assert ecog_grid[0, 0]['label'] == 'chan1'
    assert tf(trial=0, chan='chan1') == ecog_grid[0, 0]['ecog']


def simulate_data():
    chan_name = [f'chan{x + 1}' for x in range(20)]
    data = create_data(chan_name=chan_name, s_freq=1024, signal='sine', sine_freq=70, time=(0, 10))
    ampl_per_chan = array([1, 2, 5, 2, 1, 2, 3, 3, 4, 5, 1, 2, 5, 2, 1, 2, 3, 3, 4, 5])
    data.data[0] = outer(ampl_per_chan, data.data[0][0, :])  # keep same phase
    data.export(ECOG_FILE, export_format='BrainVision')
