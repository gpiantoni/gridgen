from gridgen.grid2d import make_grid_with_labels, make_grid
from gridgen.grid3d import find_vertex, construct_grid, measure_distances, measure_angles, search_grid
from gridgen.grid3d.construct import index_order
from gridgen.io import read_surf, export_grid, read_surface_ras_shift

from numpy.testing import assert_array_almost_equal
from numpy import array

from .paths import GENERATED_PATH, SMOOTH_FILE, T1_FILE


GRID_FILE = GENERATED_PATH / 'grid_020.fcsv'


def test_geometry_construct():
    offset = read_surface_ras_shift(T1_FILE)
    smooth = read_surf(SMOOTH_FILE, ras_shift=offset, normals=True)
    out_vertex = find_vertex(smooth, [-54, 3, 5])
    assert out_vertex == 37613

    vertex = search_grid(smooth, out_vertex, 5, -5)
    assert vertex == 39945

    # test construct
    grid2d = make_grid_with_labels(4, 3, 'TBLR', chan_pattern='elec{}')
    grid3d = {
        'interelec_distance': 3,
        'maximum_angle': 5,
        'step_angle': 0.2,
        }

    grid = construct_grid(smooth, out_vertex, 'elec1', grid2d['label'], grid3d, rotation=20)
    assert_array_almost_equal(
        grid['pos'][4, 2],
        array([-57.797, 1.639, -9.119]),
        decimal=3)

    offset = read_surface_ras_shift(T1_FILE)
    export_grid(grid, offset, GRID_FILE, 'slicer')
    export_grid(grid, offset, GRID_FILE, 'freeview')

    m_col, m_row = measure_distances(grid)
    assert_array_almost_equal(m_col, 3.000, decimal=3)
    assert_array_almost_equal(m_row, 3.000, decimal=3)
    m_angle = measure_angles(grid)
    assert_array_almost_equal(m_angle, 89.999, decimal=3)


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
    grid2d = make_grid(2, 2)
    grid2d['label'] = [['1', '2'], ['3', '4']]
    assert list(index_order(grid2d, '1', 'minor')) == [(0, 0), (0, 1), (1, 0), (1, 1)]
