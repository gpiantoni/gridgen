from numpy import array
from numpy.testing import assert_array_equal
from gridloc.io import read_surf, read_surface_ras_shift, export_transform, read_volume

from .paths import T1_FILE, SMOOTH_FILE, OUTPUT_PATH


def test_io_surface():
    surf = read_surf(SMOOTH_FILE, normals=False)
    assert list(surf) == ['tri', 'pos', 'tri_norm', 'pos_norm']


def test_io_volume():
    offset = read_surface_ras_shift(T1_FILE)
    assert_array_equal(offset, array([2.9809875, -6.6314774, -30.70549], dtype='f4'))
    export_transform(offset, OUTPUT_PATH / 'tkras.tfm')
