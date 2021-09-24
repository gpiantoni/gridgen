from logging import getLogger

from .utils import to_html, to_div
from .viz2d import plot_scatter, plot_grid2d
from .viz3d import plot_electrodes

lg = getLogger(__name__)

modalities = 'morphology', 'functional'


def plot_fitting(parameters, mris, model):

    output_dir = parameters['output_dir']

    scatter_file = output_dir / 'scatter'
    divs = plot_scatter(model, parameters['fit'])
    to_html(divs, scatter_file)

    grid_file = output_dir / 'projected'
    fig = plot_electrodes(mris, model['grid'], model['ecog'], functional=mris['func'])
    to_html([to_div(fig), ], grid_file)
    lg.debug(f'Exported projected to {grid_file}')

    plot_grid3d(parameters, mris, model)


def plot_grid3d(parameters, mris, model):

    output_dir = parameters['output_dir']

    fig = plot_electrodes(
        mris,
        model['grid'],
        ref_label=parameters['initial']['label'])
    grid_file = output_dir / 'electrodes.html'
    to_html([to_div(fig), ], grid_file)

    for mod in ('morphology', 'functional', 'merged'):
        grid_file = output_dir / mod
        if model.get(mod, None) is None:
            continue

        if mod == 'functional':
            func = mris['func']
        else:
            func = None

        fig1 = plot_grid2d(model[mod])
        fig2 = plot_electrodes(mris, model['grid'], model[mod], functional=func)
        if parameters.get('fit', {}).get(f'{mod}_weight', 1) < 0:
            fig1.data[-1]['reversescale'] = True
            fig2.data[-1]['marker']['reversescale'] = True

        to_html([to_div(fig1), to_div(fig2)], grid_file)
        lg.debug(f'Exported {mod} model to {grid_file}')
