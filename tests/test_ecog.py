from numpy.testing import assert_array_almost_equal

from gridloc.ecog.read_ecog import read_ecog, put_ecog_on_grid2d
from gridloc.construct import make_grid_with_labels
from gridloc.io import read_ecog2d, write_ecog2d, write_grid2d

from .paths import ECOG_FILE, OUTPUT_PATH


def test_read_ecog():
    tf = read_ecog(ECOG_FILE)
    assert tf.data[0][0] < tf.data[0][1] < tf.data[0][2]

    grid2d = make_grid_with_labels(4, 5, 'TBLR', 'chan{}')
    ecog2d = put_ecog_on_grid2d(tf, grid2d)

    assert ecog2d[0, 0]['label'] == 'chan1'
    assert tf(trial=0, chan='chan1') == ecog2d[0, 0]['ecog']

    ecog_file = OUTPUT_PATH / 'grid2d_ecog.tsv'
    grid_file = OUTPUT_PATH / 'grid2d_labels.tsv'
    write_ecog2d(ecog_file, ecog2d)
    write_grid2d(grid_file, grid2d)

    ecog2d_import = read_ecog2d(ecog_file, grid_file)
    assert_array_almost_equal(ecog2d_import['ecog'], ecog2d['ecog'])
