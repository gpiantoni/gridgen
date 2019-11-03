from scipy.optimize import basinhopping, brute

from numpy.ma import masked_invalid, corrcoef
from numpy import arange, ptp, array, sum

from .geometry import search_grid
from .morphology.distance import compute_distance
from .vascular.sphere import vascular_model
from .construct import construct_grid

def compare_models(grid, pial, angio, offset, gamma):

    m_morpho = compute_distance(grid, pial)
    m_vasc = vascular_model(grid, angio, -1 * offset)

    x = []
    for weight in arange(0.1, 1, 0.1):
        prediction = weight * normalize(m_vasc) + (1 - weight) * normalize(m_morpho)
        x.append(corrcoef_nan(-1 * prediction, gamma))

    return array(x).min()


def corrcoef_nan(A, B):
    a = masked_invalid(A)
    b = masked_invalid(B)

    msk = (~a.mask & ~b.mask)
    return corrcoef(a[msk], b[msk])[0, 1]


def normalize(x):
    return (x - x.min()) / ptp(x)

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
        niter=5,
        T=0.5,
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

    ranges = (slice(-4, 4, 4), slice(-4, 4, 4), slice(-4, 8, 4))

    res = brute(
        _compute_grid,
        ranges,
        args=args,
        disp=True,
        workers=-1,
        full_output=True,
        )

    return res
