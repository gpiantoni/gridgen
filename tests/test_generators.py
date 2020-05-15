from gridloc.construct import make_grid_with_labels


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
