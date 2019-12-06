from scipy.optimize import basinhopping, brute, minimize
from numpy.linalg import norm
from numpy.ma import masked_invalid, corrcoef
from numpy import arange, array, sum, nanmax, nanmin, argmin
from logging import getLogger

try:
    import mkl
except ImportError:
    mkl = None

from .geometry import search_grid
from .morphology.distance import compute_distance
from .vascular.sphere import vascular_model
from .construct import construct_grid
from .io import read_surf, read_surface_ras_shift

lg = getLogger(__name__)


def fitting(T1_file, dura_file, pial_file, initial, ecog, method='simplex'):

    lg.debug(f'Reading positions and computing normals of {dura_file}')
    dura = read_surf(dura_file)
    lg.debug(f'Reading positions of {pial_file}')
    pial = read_surf(pial_file, normals=False)

    offset = read_surface_ras_shift(T1_file)

    init_label = initial['label']
    init_rot = initial['rotation']
    init_ras = array(initial['RAS'])
    init_vert = argmin(norm(dura['pos'] + offset - init_ras, axis=1))
    vert_dist = norm(init_ras - offset - dura['pos'][init_vert])

    lg.debug(f'Target RAS: {init_ras}, vertex #{init_vert} RAS: {dura["pos"][init_vert]} (distance = {vert_dist:0.3}mm)')
    lg.info(f'Starting position for {init_label} is vertex #{init_vert} with orientation {initial["rotation"]}')


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
    return (x - nanmin(x)) / (nanmax(x) - nanmin(x))


def print_results(x0, res, accept):
    lg.info(f'Done: {x0[0]: 8.3f}mm {x0[1]: 8.3f}mm {x0[2]: 8.3f}° = {res: 8.3f}')


def sum_of_squares(A, B):

    a = masked_invalid(A)
    b = masked_invalid(B)

    msk = (~a.mask & ~b.mask)
    return sum((a[msk] - b[msk]) ** 2)


def _compute_grid(x0, surf, ref_vert, start_label, ecog, pial=None):
    lg.debug(f'{x0[0]: 8.3f}mm {x0[1]: 8.3f}mm {x0[2]: 8.3f}° (vert{start_vert: 6d}) = ')
    x, y, rotation = x0
    start_vert = search_grid(surf, ref_vert, x, y)
    grid = construct_grid(surf, start_vert, start_label, ecog['label'], rotation=rotation)
    model = compute_distance(grid, pial)
    SS = corrcoef_nan(ecog, model)
    lg.debug(' ' * 35 + f'{SS: 8.3f}')
    return SS


def fitting_simplex(dura, pial, ecog, start_vert, start_label, rotation):

    simplex = array([
        [-3, -3, rotation - 5],
        [3, -3, rotation - 5],
        [-3, 3, rotation - 5],
        [-3, -3, rotation + 5],
        ])

    args = (
        dura,
        start_vert,
        start_label,
        ecog,
        pial,
        )

    m = minimize(
        _compute_grid,
        array([0, 0, 0]),
        method='Nelder-Mead',
        args=args,
        options=dict(
            maxiter=100,
            initial_simplex=simplex,
            xatol=0.5,
            fatol=0.05,
            ),
        )

    return m


def fitting_hop(surf, ref_vert, start_label, grid2d, gamma, pial):
    args = (
        surf,
        ref_vert,
        start_label,
        grid2d,
        gamma,
        pial,
        )
    res = basinhopping(
        _compute_grid,
        array([0, 0, 0]),
        niter=100,
        T=0.5,
        stepsize=5,
        interval=10,
        callback=print_results,
        minimizer_kwargs=dict(
            method='Nelder-Mead',
            args=args,
            options=dict(
                maxiter=10,
                maxfev=10,
                ),
        )
    )
    return res


def fitting_brute(surf, ref_vert, start_label, grid2d, gamma, pial):
    args = (
        surf,
        ref_vert,
        start_label,
        grid2d,
        gamma,
        pial,
        )

    ranges = (
        (-10, 10),
        (-10, 10),
        (-30, 30),
        )
    if mkl is not None:
        mkl.set_num_threads(2)
    res = brute(
        _compute_grid,
        ranges,
        Ns=10,
        args=args,
        disp=True,
        workers=-1,
        full_output=True,
        finish=None,
        )

    return res
