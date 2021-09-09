"""Functions to compute the actual fitting of the grid onto the brain surface
"""
from scipy.optimize import brute, minimize
from scipy.stats import spearmanr
from multiprocessing import Pool
from numpy.linalg import norm
from numpy import arange, array, argmin, intersect1d, corrcoef
from logging import getLogger
from datetime import datetime
from json import dump

try:
    import mkl
except ImportError:
    mkl = None

from .geometry import search_grid
from .morphology.distance import compute_distance
from .vascular.sphere import compute_vasculature
from .construct import construct_grid
from .io import read_mri, write_tsv, WIRE, export_grid
from .viz import plot_results, to_div, to_html, plot_electrodes
from .examine import measure_distances, measure_angles
from .utils import be_nice, match_labels, normalize

lg = getLogger(__name__)


def fitting(output, ecog, grid3d, initial, mri, fit):
    """Fit the brain activity onto the surface

    Parameters
    ----------
    output : path
        folder to export to
    ecog : numpy 2d array
        array with labels and ecog values
    grid3d : dict
        see parameters.md
    initial : dict
        see parameters.md
    mri : dict
        see parameters.md
    fit : dict
        see parameters.md

    Returns
    -------
    instance of grid2d
        grid2d with best positions
    """
    start_time = datetime.now()
    mris = read_mri(**mri)

    init_ras = array(initial['RAS'])
    initial['vertex'] = argmin(norm(mris['dura']['pos'] - init_ras, axis=1))
    vert_dist = norm(init_ras - mris['dura']['pos'][initial["vertex"]])

    lg.debug(f'Target RAS: {init_ras}, vertex #{initial["vertex"]} RAS: {mris["dura"]["pos"][initial["vertex"]]} (distance = {vert_dist:0.3}mm)')
    lg.info(f'Starting position for {initial["label"]} is vertex #{initial["vertex"]} with orientation {initial["rotation"]}')

    # has to be a tuple
    minimizer_args = (
        ecog,  # 0
        grid3d,  # 1
        initial,  # 2
        mris,  # 3
        fit,  # 4
        )

    if not fit:
        # start position
        grid = construct_grid(
            mris["dura"],
            initial['vertex'],
            initial["label"],
            ecog['label'],
            grid3d,
            rotation=initial['rotation'])
        fig = plot_electrodes(mris['pial'], grid, ref_label=initial['label'], angio=mris['angio'])
        grid_file = output / 'start_pos.html'
        to_html([to_div(fig), ], grid_file)
        return

    if fit['method'] == 'simplex':
        m = fitting_simplex(corr_ecog_model, None, minimizer_args)
        best_fit = m.x

    elif fit['method'] == 'brute':
        m = fitting_brute(corr_ecog_model, minimizer_args)
        best_fit = m[0]

    end_time = datetime.now()
    comp_dur = (end_time - start_time).total_seconds()
    lg.debug(f'Model fitting took {comp_dur:1.0f}s')

    # create grid with best values
    x, y, rotation = best_fit
    model = corr_ecog_model(best_fit, *minimizer_args, final=True)
    lg.info(f'Best fit at {x:+8.3f}mm {y:+8.3f}mm {rotation:+8.3f}° (vert{model["vert"]: 6d}) = {model["cc"]:+8.3f} (# included channels:{len(model["n_chan"]): 4d}, vascular contribution: {model["percent_vasc"]:.2f}%)')

    plot_results(model, mris['pial'], mris['ras_shift'], output, angio=mris['angio'])

    model = remove_wires(model)

    measure_distances(model['grid'])
    measure_angles(model['grid'])
    out = {
        'label': initial['label'],
        'surface': str(mri['dura_file']),
        'vertex': int(model['vert']),
        'pos': list(mris['dura']['pos'][model['vert'], :]),
        'normals': list(mris['dura']['pos_norm'][model['vert'], :]),
        'rotation': rotation,
        'percent_vasculature': model['percent_vasc'],
        'n_included_channels': model['n_chan'],
        'corrcoef': model['cc'],
        'duration': comp_dur,
        }
    results_file = output / 'results.json'
    with results_file.open('w') as f:
        dump(out, f, indent=2)

    grid_file = output / 'electrodes'
    export_grid(model['grid'], mris['ras_shift'], grid_file)
    write_tsv(model['grid']['label'], model['grid']['pos'], grid_file)
    lg.debug(f'Exported electrodes to {grid_file} (coordinates in MRI volume space, not mesh space)')

    return model['grid']


def corr_ecog_model(x0, ecog, grid3d, initial, mri, fit, final=False):
    """Main model to minimize

    Note that when final=False, cc should be minimized (the smaller the better)
    while when final=True, cc should be as large as possible (more intuitive
    reading)

    Parameters
    ----------
    ranges : None
        ignored, but easier to pass when working with arguments

    """
    x, y, rotation = x0
    start_vert = search_grid(mri["dura"], initial["vertex"], x, y)
    grid = construct_grid(
        mri["dura"],
        start_vert,
        initial["label"],
        ecog['label'],
        grid3d,
        rotation=initial['rotation'] + rotation)

    morpho = compute_distance(grid, mri['pial'], fit['morphology'], fit['maximum_distance'])

    if mri['angio'] is not None:
        vasc = compute_vasculature(grid, mri['angio'])
        e, m, v = match_labels(ecog, morpho, vasc)

    else:
        e, m = match_labels(ecog, morpho)
        v = vasc = None

    i, cc = compare_models(e, m, v, correlation=fit['correlation'])
    if not final:
        lg.debug(f'{x0[0]:+8.3f}mm {x0[1]:+8.3f}mm {x0[2]:+8.3f}° (vert{start_vert: 6d}) = {cc * -1:+8.3f} (# included channels:{len(e): 4d}, vascular contribution: {100 * (1 - i):.2f}%)')

    if final:
        return {
            'ecog': ecog,
            'vert': start_vert,
            'grid': grid,
            'morpho': morpho,
            'vasc': vasc,
            'percent_vasc': 100 * (1 - i),
            'n_chan': len(e),  # number of channels used to compute
            'cc': -1 * cc,
            }

    else:
        return cc


def compare_models(E, M, V=None, correlation='parametric'):

    E = normalize(E)
    M = normalize(M)

    if V is not None:
        V = normalize(V)
        WEIGHTS = arange(0, 1.1, 0.1)
    else:
        V = 0
        WEIGHTS = [1, ]

    x = []
    for weight in WEIGHTS:
        prediction = weight * M + (1 - weight) * V

        if correlation == 'parametric':
            c = corrcoef(E, prediction)[0, 1]
        else:
            c = spearmanr(E, prediction).correlation

        x.append(c)

    x = array(x)
    i = argmin(x)

    return WEIGHTS[i], x[i]


def corrcoef_match(ecog, estimate, field='morphology'):
    """correlation but make sure that the labels match
    """
    good = ecog['label'][ecog['good']]
    ecog_id = intersect1d(ecog['label'], good, return_indices=True)[1]
    a = ecog['ecog'].flatten('C')[ecog_id]

    estimate_id = intersect1d(estimate['label'], good, return_indices=True)[1]
    b = estimate[field].flatten('C')[estimate_id]

    return corrcoef(a, b)[0, 1]


def fitting_brute(func, args):

    ranges = args[4]['ranges']
    # make sure that the last point is included in the range
    for k, v in ranges.items():
        ranges[k][1], ranges[k][2] = ranges[k][2], ranges[k][1]
        ranges[k][1] += ranges[k][2]

    ranges = (
        slice(ranges['x']),
        slice(ranges['y']),
        slice(ranges['rotation']),
        )

    if mkl is not None:
        mkl.set_num_threads(2)

    with Pool(initializer=be_nice) as p:
        res = brute(
            corr_ecog_model,
            ranges,
            args=args,
            disp=True,
            workers=p.map,
            full_output=True,
            finish=fitting_simplex,
            )

    return res


def fitting_simplex(func, init, args):

    if init is None:  # when called stand alone
        x = y = rotation = 0
        steps = args[4]['steps']
    else:
        lg.info(f'Applying simplex from starting point: {init[0]:+8.3f}mm {init[1]:+8.3f}mm {init[2]:+8.3f}°')
        x, y, rotation = init
        # convert ranges to simplex steps
        steps = {k: args[4]['ranges'][k][2] / 2 for k in ('x', 'y', 'rotation')}

    simplex = array([
        [x - steps['x'], y - steps['y'], rotation - steps['rotation']],
        [x + steps['x'], y - steps['y'], rotation - steps['rotation']],
        [x - steps['x'], y + steps['y'], rotation - steps['rotation']],
        [x - steps['x'], y - steps['y'], rotation + steps['rotation']],
        ])

    m = minimize(
        func,
        array([0, 0, 0]),  # ignored
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


def remove_wires(model):
    """We select rows with WIRE, so that we keep the 2d shape of the models
    """
    i_keep = model['grid']['label'] != WIRE
    i_keep = i_keep.all(axis=1)
    for field in ['ecog', 'grid', 'morpho', 'vasc']:
        if model[field] is not None:
            model[field] = model[field][i_keep, :]

    return model
