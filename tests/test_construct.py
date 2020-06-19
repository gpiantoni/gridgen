from gridloc.construct import construct_grid, make_grid_with_labels
from gridloc.io import read_surf, export_grid
from gridloc.geometry import search_grid

from numpy.testing import assert_array_almost_equal
from numpy import array

from .paths import GENERATED_PATH, SMOOTH_FILE


def test_geometry_construct():
    # test geometry
    smooth = read_surf(SMOOTH_FILE, normals=True)
    out_vertex = search_grid(smooth, 30000, 5, 5)
    assert out_vertex == 27729

    # test construct
    grid_file = GENERATED_PATH / 'grid_020.fcsv'

    grid2d = make_grid_with_labels(8, 8, chan_pattern='elec{}')
    grid = construct_grid(smooth, 33154, 'elec1', grid2d['label'], rotation=20)
    assert_array_almost_equal(
        grid['pos'][0, 0],
        array([7.779, 0.410, 52.968]),
        decimal=3)

    export_grid(grid, grid_file, 'slicer')
    export_grid(grid, grid_file, 'freeview')
