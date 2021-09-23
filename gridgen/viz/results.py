from logging import getLogger

from .utils import to_html, to_div
from .viz2d import plot_scatter, plot_grid2d
from .viz3d import plot_electrodes

lg = getLogger(__name__)

modalities = 'morphology', 'functional'


def plot_grid3d(parameters, mris, model):
    output_dir = parameters['output_dir']

    fig = plot_electrodes(
        mris,
        model['grid'],
        ref_label=parameters['initial']['label'])
    grid_file = output_dir / 'electrodes.html'
    to_html([to_div(fig), ], grid_file)

    _plot_2d_and_3d(output_dir, mris, model)


def plot_fitting(output_dir, mris, model, fit):

    scatter_file = output_dir / 'scatter'
    divs = plot_scatter(model, fit)
    to_html(divs, scatter_file)

    grid_file = output_dir / 'projected'
    fig = plot_electrodes(mris, model['grid'], model['ecog'], functional=mris['func'])
    to_html([to_div(fig), ], grid_file)
    lg.debug(f'Exported projected to {grid_file}')

    _plot_2d_and_3d(output_dir, mris, model, fit)


def _plot_2d_and_3d(output_dir, mris, model, fit={}):
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
        if fit.get(f'{mod}_weight', '') == 'negative':
            fig1.data[-1]['reversescale'] = True
            fig2.data[-1]['marker']['reversescale'] = True

        to_html([to_div(fig1), to_div(fig2)], grid_file)
        lg.debug(f'Exported {mod} model to {grid_file}')
