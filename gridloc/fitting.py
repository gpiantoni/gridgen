from scipy.optimize import basinhopping, brute
from numpy import array, sum
from .geometry import search_grid
from .morphology.distance import compute_distance
from .construct import construct_grid


def print_results(x0, res, accept):
    print(x0, res)


def sum_of_squares(A, B):
    return sum((A - B) ** 2)


def _compute_grid(x0, surf, ref_vert, start_label, n_rows, n_cols, gamma,
                  pial=None):

    print('.', end='')
    x, y, rotation = x0
    start_vert = search_grid(surf, 31859, x, y)
    grid = construct_grid(surf, start_vert, start_label, n_rows, n_cols, rotation=rotation)
    model = compute_distance(grid, pial)
    return sum_of_squares(gamma, model)


def fitting_hop(surf, ref_vert, start_label, n_rows, n_cols, gamma, pial):
    args = (
        surf,
        ref_vert,
        start_label,
        n_rows,
        n_cols,
        gamma,
        pial,
        )
    res = basinhopping(
        _compute_grid,
        array([1, 1, 3]),
        niter=100,
        T=1.0,
        stepsize=1,
        callback=print_results,
        minimizer_kwargs=dict(
            method='Nelder-Mead',
            args=args,
            options=dict(
                maxiter=10,
                ),
        )
    )
    return res


def fitting_brute(surf, ref_vert, start_label, n_rows, n_cols, gamma, pial):
    args = (
        surf,
        ref_vert,
        start_label,
        n_rows,
        n_cols,
        gamma,
        pial,
        )

    ranges = (slice(-4, 4, 2), slice(-4, 4, 2), slice(-4, 8, 4))

    res = brute(
        _compute_grid,
        ranges,
        args=args,
        disp=True,
        workers=-1,
        full_output=True,
        )

    return res
