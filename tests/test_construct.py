from gridloc.construct import construct_grid, make_grid_with_labels, index_order
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


def test_make_grids_with_labels():
    n_rows = 4
    n_columns = 3
    chan_pattern = 'chan{:02d}'

    grid2d = make_grid_with_labels(n_rows, n_columns, 'TBLR', chan_pattern)
    assert grid2d['label'][0, 0] == chan_pattern.format(1)  # TL
    assert grid2d['label'][n_rows - 1, 0] == chan_pattern.format(n_rows)  # TR
    assert grid2d['label'][0, n_columns - 1] == chan_pattern.format((n_rows - 1) * n_columns)  # BL
    assert grid2d['label'][n_rows - 1, n_columns - 1] == chan_pattern.format(n_rows * n_columns)  # BR

    grid2d = make_grid_with_labels(n_rows, n_columns, 'RLBT', chan_pattern)
    assert grid2d['label'][0, 0] == chan_pattern.format(n_rows * n_columns)  # TL
    assert grid2d['label'][n_rows - 1, 0] == chan_pattern.format(n_columns)  # TR
    assert grid2d['label'][0, n_columns - 1] == chan_pattern.format((n_rows - 1) * n_columns + 1)  # BL
    assert grid2d['label'][n_rows - 1, n_columns - 1] == chan_pattern.format(1)  # BR


def test_index_order():
    grid2d = make_grid_with_labels(2, 2)
    assert list(index_order(grid2d, '1', 'minor')) == [(0, 0), (0, 1), (1, 0), (1, 1)]
