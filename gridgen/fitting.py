"""Functions to compute the actual fitting of the grid onto the brain surface
"""
from scipy.optimize import brute, minimize
from scipy.stats import spearmanr
from multiprocessing import Pool
from numpy import arange, array, argmin, corrcoef, zeros, dtype, intersect1d, NaN, unravel_index
from logging import getLogger
from datetime import datetime
from json import dump

try:
    import mkl
except ImportError:
    mkl = None

from .grid3d import construct_grid
"""
from .geometry import search_grid
from .morphology.distance import compute_distance
from .vascular.sphere import compute_vasculature
from .io import read_mri, write_tsv, WIRE, export_grid
from .viz import plot_results, to_div, to_html, plot_electrodes
from .examine import measure_distances, measure_angles
from .utils import be_nice, match_labels, normalize
"""

lg = getLogger(__name__)


def make_grid3d_model(output, grid2d, mris, grid3d, initial, morphology={}, functional={}):

    grid = construct_grid(
        mris["dura"],
        initial['vertex'],
        initial["label"],
        grid2d['label'],
        grid3d,
        rotation=initial['rotation'])

    model = compute_model(mris, grid, morphology, functional)

    fig = plot_electrodes(mris['pial'], grid, ref_label=initial['label'], angio=mris['angio'])
    grid_file = output / 'start_pos.html'
    to_html([to_div(fig), ], grid_file)



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

    lg.info(f'Starting position for {initial["label"]} is vertex #{initial["vertex"]} with orientation {initial["rotation"]}')

    # has to be a tuple
    minimizer_args = (
        ecog,  # 0
        grid3d,  # 1
        initial,  # 2
        mris,  # 3
        fit,  # 4
        )

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
    lg.info(f'Best fit at {x:+8.3f}mm {y:+8.3f}mm {rotation:+8.3f}° (vert{model["vertex"]: 6d}) = {model["corr_coef"]:+8.3f} (# included channels:{model["n_channels"]: 4d}, vascular contribution: {model["percent_vasc"]:.2f}%)')

    plot_results(model, mris['pial'], output, angio=mris['angio'])

    model = remove_wires(model)

    measure_distances(model['grid'])
    measure_angles(model['grid'])
    out = {
        'label': initial['label'],
        'surface': str(mri['dura_file']),
        'vertex': int(model['vertex']),
        'pos': list(mris['dura']['pos'][model['vertex'], :]),
        'normals': list(mris['dura']['pos_norm'][model['vertex'], :]),
        'rotation': rotation,
        'percent_vasculature': model['percent_vasc'],
        'n_included_channels': model['n_channels'],
        'corr_coef': model['corr_coef'],
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

    morpho = compute_distance(grid, mri['pial'], fit['distance'], fit['maximum_distance'])
    morpho['value'] = morpho['value'] ** (-1 * fit['penalty'])

    if mri['angio'] is not None:
        vasc = compute_vasculature(grid, mri['angio'])
        e, m, v = match_labels(ecog, morpho, vasc)[1:]

    else:
        e, m = match_labels(ecog, morpho)[1:]
        v = vasc = None

    i, cc = compare_models(e, m, v, correlation=fit['correlation'])
    if not final:
        lg.debug(f'{x0[0]:+8.3f}mm {x0[1]:+8.3f}mm {x0[2]:+8.3f}° (vert{start_vert: 6d}) = {cc:+8.3f} (# included channels:{len(e): 4d}, vascular contribution: {100 * (1 - i):.2f}%)')

    if final:
        model = {
            'ecog': ecog,
            'vertex': start_vert,
            'grid': grid,
            'morphology': morpho,
            'vasculature': vasc,
            'percent_vasc': 100 * (1 - i),
            'n_channels': len(e),  # number of channels used to compute
            'corr_coef': cc,
            }
        model['merged'] = merge_models(model)
        return model

    else:
        return -1 * cc  # value to minimize


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


def merge_models(model):

    d_ = dtype([
        ('label', '<U256'),   # labels cannot be longer than 256 char
        ('value', 'f4'),
        ])
    merged = zeros((model['ecog'].shape[0], model['ecog'].shape[1]), dtype=d_)
    merged['label'] = model['ecog']['label']

    if model['vasculature'] is None:
        merged['value'] = normalize(model['ecog']['value'])

    else:
        merged['value'].fill(NaN)

        labels, e, m, v = match_labels(
            model['ecog'],
            model['morphology'],
            model['vasculature']
            )
        percent_vasc = model['percent_vasc']
        M = normalize(m)
        V = normalize(v)
        prediction = V * percent_vasc / 100 + M * (100 - percent_vasc) / 100

        [i0, i1] = intersect1d(merged['label'], labels, return_indices=True)[1:]
        i0r, i0c = unravel_index(i0, merged.shape)
        merged['value'][i0r, i0c] = prediction[i1]

    return merged
