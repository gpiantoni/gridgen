"""Functions to compute the actual fitting of the grid onto the brain surface
"""
from scipy.optimize import brute, minimize
from multiprocessing import Pool
from numpy import array, mean
from logging import getLogger
from datetime import datetime
from json import dump

try:
    import mkl
except ImportError:
    mkl = None

from .io import export_electrodes
from .grid3d import construct_grid, search_grid, find_vertex, measure_distances, measure_angles
from .models import compute_model, merge_models, compare_model_with_ecog, sum_models
from .utils import be_nice, remove_wires, _JSONEncoder_path
from .viz import plot_fitting, plot_grid3d

lg = getLogger(__name__)
MAXITER = 200


def fitting(output, ecog, mris, grid3d, initial, fit, morphology, functional):
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
    morphology : dict
        see parameters.md
    functional : dict
        see parameters.md

    Returns
    -------
    instance of grid2d
        grid2d with best positions
    """
    start_time = datetime.now()
    initial["vertex"] = find_vertex(mris['dura'], initial['RAS'])
    lg.info(f'Starting position for {initial["label"]} is vertex #{initial["vertex"]} with orientation {initial["rotation"]}')

    params = {
        'initial': initial,
        'grid3d': grid3d,
        'morphology': morphology,
        'functional': functional,
        'fit': fit,
        }

    # has to be a tuple
    minimizer_args = (
        ecog,  # 0
        mris,  # 1
        params,  # 2
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
    lg.info(f'Best fit at {x:+8.3f}mm {y:+8.3f}mm {rotation:+8.3f}° (vert{model["vertex"]: 6d}) = {model["summary"]:+8.3f} (# included channels:{model["n_channels"]: 4d}, functional contribution: {model["percent_functional"]:.2f}%)')

    params['output_dir'] = output

    if fit['metric'] == 'sum':
        plot_grid3d(params, mris, model)
    else:
        plot_fitting(params, mris, model)

    model = remove_wires(model)

    out = {
        'label': initial['label'],
        'vertex': model['vertex'],
        'pos': list(mris['dura']['pos'][model['vertex'], :]),
        'normals': list(mris['dura']['pos_norm'][model['vertex'], :]),
        'rotation': initial['rotation'] + rotation,
        'percent_functional': model['percent_functional'],
        'n_included_channels': model['n_channels'],
        'summary': model['summary'],
        'duration': comp_dur,
        'mean_elec_distance': mean(measure_distances(model['grid'])),
        'mean_angle': measure_angles(model['grid']),
        }
    results_file = output / 'results.json'
    with results_file.open('w') as f:
        dump(out, f, indent=2, cls=_JSONEncoder_path)

    export_electrodes(output, model, mris)

    return model


def corr_ecog_model(x0, ecog, mris, params, final=False):
    """Main model to minimize


    Parameters
    ----------
    x0 : list of 3 floats
        start point to look for vertex
    ecog : (n_rows, n_cols) array
        array with ecog values (can be None when fit['metric'] is sum)
    mris : dict
        MRIs and meshes useful for computing the model
    params : dict
        initial, grid3d, morphology, functional
    final : bool
        whether to compute only the correlation value (False) or the full model
        (True)

    Returns
    -------
    float (when final=False)
        value to minimize
    dict (when final=True)
        full model

    Note
    ----
    When final=False, cc should be minimized (the smaller the better)
    while when final=True, cc should be as large as possible (more intuitive
    reading)
    """
    x, y, rotation = x0
    start_vertex = search_grid(mris["dura"], params['initial']["vertex"], x, y)

    grid = construct_grid(
        mris["dura"],
        start_vertex,
        params['initial']["label"],
        ecog['label'],
        params['grid3d'],
        rotation=params['initial']['rotation'] + rotation)

    model = compute_model(mris, grid, params['morphology'], params['functional'])

    if params['fit']['metric'] == 'sum':
        weight, cc, chans, vals = sum_models(model, fit=params['fit'])

    else:
        weight, cc, chans, vals = compare_model_with_ecog(model, ecog, fit=params['fit'])

    if not final:
        lg.debug(f'{x0[0]:+8.3f}mm {x0[1]:+8.3f}mm {x0[2]:+8.3f}° (vert{start_vertex: 6d}) = {cc:+8.3f} (# included channels:{len(chans): 4d}, functional contribution: {weight:.2f}%)')
        return -1 * cc  # value to minimize

    else:
        model['ecog'] = ecog
        model['vertex'] = start_vertex
        model['percent_functional'] = weight
        model['summary'] = cc
        model['n_channels'] = len(chans)
        model['merged'] = merge_models(ecog, chans, vals)
        return model


def fitting_brute(func, args):
    """Run the fitting procedure using the brute approach

    Parameters
    ----------
    func : func
        function to minimize
    args : tuple
        rest of the parameters

    Returns
    -------
    output of the model
    """
    ranges = args[2]['fit']['ranges']
    # make sure that the last point is included in the range
    for k, v in ranges.items():
        ranges[k][1], ranges[k][2] = ranges[k][2], ranges[k][1]
        ranges[k][1] += ranges[k][2]

    ranges = (
        slice(*ranges['x']),
        slice(*ranges['y']),
        slice(*ranges['rotation']),
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
    """Run the fitting procedure using the simplex approach

    Parameters
    ----------
    func : func
        function to minimize
    args : tuple
        rest of the parameters

    Returns
    -------
    output of the model
    """

    if init is None:  # when called stand alone
        x = y = rotation = 0
        steps = args[2]['fit']['steps']
    else:
        lg.info(f'Applying simplex from starting point: {init[0]:+8.3f}mm {init[1]:+8.3f}mm {init[2]:+8.3f}°')
        x, y, rotation = init
        # convert ranges to simplex steps
        steps = {k: args[2]['fit']['ranges'][k][2] / 2 for k in ('x', 'y', 'rotation')}

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
            maxiter=MAXITER,
            initial_simplex=simplex,
            xatol=0.5,
            fatol=0.05,
            ),
        )

    return m
