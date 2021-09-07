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
from .io import read_surf, read_surface_ras_shift, read_volume, write_tsv
from .viz import plot_results, to_div, to_html, plot_electrodes
from .examine import measure_distances, measure_angles
from .utils import be_nice, match_labels, normalize

lg = getLogger(__name__)


def fitting(T1_file, dura_file, pial_file, initial, ecog, output, angio_file=None,
            angio_threshold=None, correlation='parametric', ranges={}, steps={},
            method='brute'):
    """Fit the brain activity onto the surface

    Parameters
    ----------
    T1_file : path
        path to T1 image (in particular, the T1.mgz from freesurfer)
    dura_file : path
        path to dura surface (for example, the smoothed pial surface)
    pial_file : path
        path to pial surface (in particular, the lh.pial or rh.pial from freesurfer)
    initial : dict
        start position for search, with fields:
            - label : label for the reference electrode
            - RAS : initial location for the reference electrode
            - rotation : degree of rotation of the grid (in degrees, 0° is roughly pointing up)
    angio_file : path or None
        path to angiogram (in NIfTI format). Optional.
    angio_threshold : float
        value to threshold the angio_file. Optional.
    correlation : str
        'parametric' (Pearson) or 'nonparametric' (rank)
    method : str
        'simplex', 'brute'
    ranges : dict of lists
        keys are x-direction, y-direction, rotation
    steps : dict of lists
        keys are x-direction, y-direction, rotation

    Returns
    -------
    instance of grid2d
        grid2d with best positions
    """
    start_time = datetime.now()

    ras_shift = read_surface_ras_shift(T1_file)
    lg.debug(f'Reading positions and computing normals of {dura_file}')
    dura = read_surf(dura_file, ras_shift=ras_shift)
    lg.debug(f'Reading positions of {pial_file}')
    pial = read_surf(pial_file, normals=False, ras_shift=ras_shift)

    if angio_file is not None and angio_file:
        lg.debug(f'Reading angiogram from {angio_file} and thresholding at {angio_threshold}')
        angio = read_volume(angio_file, angio_threshold)
    else:
        angio = None

    ref_label = initial['label']
    init_rot = initial['rotation']
    init_ras = array(initial['RAS'])
    init_vert = argmin(norm(dura['pos'] - init_ras, axis=1))
    vert_dist = norm(init_ras - dura['pos'][init_vert])

    lg.debug(f'Target RAS: {init_ras}, vertex #{init_vert} RAS: {dura["pos"][init_vert]} (distance = {vert_dist:0.3}mm)')
    lg.info(f'Starting position for {ref_label} is vertex #{init_vert} with orientation {initial["rotation"]}')

    # start position is init_vert, plus some rotation
    x = y = rotation = 0

    # make sure that the last point is included in the range
    for k, v in ranges.items():
        ranges[k][1], ranges[k][2] = ranges[k][2], ranges[k][1]
        ranges[k][1] += ranges[k][2]

    # has to be a tuple
    minimizer_args = (
        dura,  # 0
        init_vert,  # 1
        ref_label,  # 2
        init_rot,  # 3
        ecog,  # 4
        pial,  # 5
        angio,  # 6
        correlation,  # 7
        ranges,  # 8
        )

    # start position
    model = corr_ecog_model([x, y, rotation], *minimizer_args[:-1], final=True)
    lg.info(f'Start position at {x: 8.3f}mm {y: 8.3f}mm {rotation: 8.3f}° (vert{model["vert"]: 6d}) = {model["cc"]: 8.3f} (vascular contribution: {model["percent_vasc"]:.2f}%)')
    fig = plot_electrodes(pial, model['grid'], ref_label=ref_label)
    grid_file = output / 'start_pos.html'
    to_html([to_div(fig), ], grid_file)

    if method == 'simplex':
        m = fitting_simplex(corr_ecog_model, [x, y, rotation], minimizer_args)
        best_fit = m.x

    elif method == 'brute':
        m = fitting_brute(corr_ecog_model, [x, y, rotation], minimizer_args)
        best_fit = m[0]

    end_time = datetime.now()
    comp_dur = (end_time - start_time).total_seconds()
    lg.debug(f'Model fitting took {comp_dur:1.0f}s')

    # create grid with best values
    x, y, rotation = best_fit
    model = corr_ecog_model(best_fit, *minimizer_args[:-1], final=True)
    lg.info(f'Best fit at {x: 8.3f}mm {y: 8.3f}mm {rotation: 8.3f}° (vert{model["vert"]: 6d}) = {model["cc"]: 8.3f} (vascular contribution: {model["percent_vasc"]:.2f}%)')

    measure_distances(model['grid'])
    measure_angles(model['grid'])

    out = {
        'ref_label': ref_label,
        'surface': str(dura_file),
        'vertex': int(model['vert']),
        'pos': list(dura['pos'][model['vert'], :]),
        'normals': list(dura['pos_norm'][model['vert'], :]),
        'rotation': rotation,
        'percent_vasculature': model['percent_vasc'],
        'r2': model['cc'],
        'duration': comp_dur,
        }
    results_file = output / 'results.json'
    with results_file.open('w') as f:
        dump(out, f, indent=2)

    grid_file = output / 'electrodes'
    write_tsv(model['grid']['label'], model['grid']['pos'], grid_file)
    lg.debug(f'Exported electrodes to {grid_file} (coordinates in MRI volume space, not mesh space)')

    plot_results(model, pial, ras_shift, output)
    return model['grid']


def corr_ecog_model(x0, dura, ref_vert, ref_label, init_rot, ecog, pial, angio=None,
                    correlation=None, ranges=None, final=False):
    """Main model to minimize

    Parameters
    ----------
    ranges : None
        ignored, but easier to pass when working with arguments

    """
    x, y, rotation = x0
    start_vert = search_grid(dura, ref_vert, x, y)
    grid = construct_grid(dura, start_vert, ref_label, ecog['label'], rotation=init_rot + rotation)

    morpho = compute_distance(grid, pial, 'minimum')

    if angio is not None and angio is not False:
        vasc = compute_vasculature(grid, angio)
        e, m, v = match_labels(ecog, morpho, vasc)

    else:
        e, m = match_labels(ecog, morpho)
        v = vasc = None

    i, cc = compare_models(e, m, v, correlation=correlation)
    if not final:
        lg.debug(f'{x0[0]: 8.3f}mm {x0[1]: 8.3f}mm {x0[2]: 8.3f}° (vert{start_vert: 6d}) = {cc: 8.3f} (vascular contribution: {100 * (1 - i):.2f}%)')

    if final:
        return {
            'ecog': ecog,
            'vert': start_vert,
            'grid': grid,
            'morpho': morpho,
            'vasc': vasc,
            'percent_vasc': 100 * (1 - i),
            'cc': cc,
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


def fitting_brute(func, init, args):

    # shift rotation by initial value
    # for x and y, the default value MUST be zero, because we start at
    # the reference vertex
    rotation = args[8]['rotation']
    rotation[0] += init[2]
    rotation[1] += init[2]

    ranges = (
        slice(*args[7]['x']),
        slice(*args[7]['y']),
        slice(*rotation),
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
    """TODO: check after brute"""
    lg.info(f'Starting point: {init[0]: 8.3f}mm {init[1]: 8.3f}mm {init[2]: 8.3f}°. Now applying simplex')

    x, y, rotation = init
    steps = args[7]
    simplex = array([
        [x - steps['x'], y - steps['y'], rotation - ranges['rotation'][2] / 2],
        [x + steps['x'], y - steps['y'], rotation - ranges['rotation'][2] / 2],
        [x - steps['x'], y + steps['y'], rotation - ranges['rotation'][2] / 2],
        [x - steps['x'], y - steps['y'], rotation + ranges['rotation'][2] / 2],
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
