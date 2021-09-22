from numpy import array
from numpy.testing import assert_array_almost_equal

from gridgen.io import read_surf, read_surface_ras_shift, export_transform, read_volume

from .paths import T1_FILE, SMOOTH_FILE, OUTPUT_PATH


def test_io_surface():
    surf = read_surf(SMOOTH_FILE, normals=False)
    assert list(surf) == ['tri', 'pos', 'tri_norm', 'pos_norm']


def test_io_volume():
    offset = read_surface_ras_shift(T1_FILE)
    assert_array_almost_equal(
        offset,
        array([3.663, 12.814, -31.671], dtype='f4'),
        decimal=3,
        )
    export_transform(offset, OUTPUT_PATH / 'tkras.tfm')

    T1 = read_volume(T1_FILE, 150)
    assert T1.shape[0] == 862
