from logging import getLogger

from .utils import to_html, to_div
from .viz2d import plot_scatter, plot_grid2d
from .viz3d import plot_electrodes

lg = getLogger(__name__)

modalities = 'morphology', 'functional'


def plot_grid3d(parameters, mris, model):
    output_dir = parameters['output_dir']

    fig = plot_electrodes(
        mris['pial'],
        model['grid'],
        ref_label=parameters['initial']['label'])
    grid_file = output_dir / 'electrodes.html'
    to_html([to_div(fig), ], grid_file)

    for mod in modalities:
        grid_file = output_dir / mod
        if mod == 'functional':
            func = mris['func']
        else:
            func = None
        fig1 = plot_grid2d(model[mod])
        fig2 = plot_electrodes(mris['pial'], model['grid'], model[mod], functional=func)
        to_html([to_div(fig1), to_div(fig2)], grid_file)
        lg.debug(f'Exported {mod} model to {grid_file}')


def plot_fitting(parameters, mris, model):

    scatter_file = output / 'scatter'
    divs = plot_scatter(model)
    to_html(divs, scatter_file)

    grid_file = output / 'projected'
    fig = plot_electrodes(pial, model, 'ecog', angio=angio)
    to_html([to_div(fig), ], grid_file)
    lg.debug(f'Exported merged model to {grid_file}')

    grid_file = output / 'morphology'
    fig0 = plot_grid2d(model['morphology'], 'morphology')
    fig1 = plot_electrodes(pial, model, 'morphology')
    to_html([to_div(fig0), to_div(fig1)], grid_file)

    if model['functional'] is not None:
        grid_file = output / 'functional'
        fig0 = plot_grid2d(model['functional'], 'functional')
        fig1 = plot_electrodes(pial, model, 'functional', angio=angio)
        to_html([to_div(fig0), to_div(fig1)], grid_file)
        lg.debug(f'Exported vascular to {grid_file}')

        grid_file = output / 'merged'
        fig0 = plot_grid2d(model['merged'], 'merged')
        fig1 = plot_electrodes(pial, model, 'merged', angio=angio)
        to_html([to_div(fig0), to_div(fig1)], grid_file)
        lg.debug(f'Exported merged model to {grid_file}')
