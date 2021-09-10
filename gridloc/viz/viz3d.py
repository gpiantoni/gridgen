from logging import getLogger
from numpy import mean, sign, nanmin, nanmax, log10
import plotly.graph_objects as go

from ..io import WIRE
from .utils import to_html, to_div, default_colorbar, COLORSCALE
from .viz2d import plot_scatter, plot_grid2d


AXIS = dict(
    title="",
    visible=False,
    zeroline=False,
    showline=False,
    showticklabels=False,
    showgrid=False,
    )
MARKER_SIZE = 3


lg = getLogger(__name__)


def plot_results(model, pial, output, angio=None):

    scatter_file = output / 'scatter'
    divs = plot_scatter(model)
    to_html(divs, scatter_file)

    grid_file = output / 'projected'
    fig = plot_electrodes(pial, model, 'ecog', angio=angio)
    to_html([to_div(fig), ], grid_file)
    lg.debug(f'Exported merged model to {grid_file}')

    grid_file = output / 'morphology'
    fig0 = plot_grid2d(model['morpho'], 'morphology')
    fig1 = plot_electrodes(pial, model, 'morphology')
    to_html([to_div(fig0), to_div(fig1)], grid_file)

    if model['vasc'] is not None:
        grid_file = output / 'vascular'
        fig0 = plot_grid2d(model['vasc'], 'vasculature')
        fig1 = plot_electrodes(pial, model, 'vasculature', angio=angio)
        to_html([to_div(fig0), to_div(fig1)], grid_file)
        lg.debug(f'Exported vascular to {grid_file}')

        grid_file = output / 'merged'
        fig0 = plot_grid2d(model['merged'], 'merged')
        fig1 = plot_electrodes(pial, model, 'merged', angio=angio)
        to_html([to_div(fig0), to_div(fig1)], grid_file)
        lg.debug(f'Exported merged model to {grid_file}')


def plot_electrodes(pial, model, value=None, ref_label=None, angio=None):
    grid = model['grid']
    values = model[value]
    right_or_left = sign(mean(pial['pos'][:, 0]))
    pos = grid['pos'].reshape(-1, 3)
    norm = grid['norm'].reshape(-1, 3)
    labels = grid['label'].reshape(-1)

    if values is None:
        iswire = labels == WIRE
        colors = labels.copy()
        colors[iswire] = 'red'
        colors[~iswire] = 'black'
        if ref_label is not None:
            colors[labels == ref_label] = 'green'
        marker = dict(
            size=MARKER_SIZE,
            color=colors,
            )
        hovertext = labels

    else:

        colorbar = default_colorbar(value)
        values = values.reshape(-1)
        marker = dict(
            size=MARKER_SIZE,
            color=values,
            colorscale=COLORSCALE,
            showscale=True,
            cmin=nanmin(values),
            cmax=nanmax(values),
            colorbar=dict(
                title=dict(
                    text=colorbar.replace(' ', '<br>'),
                    )
                ),
            )
        hovertext = [f'{x0}<br>{x1:0.3f}' for x0, x1 in zip(labels, values)]

    if value == 'morphology':
        traces = [
            go.Cone(
                x=pos[:, 0],
                y=pos[:, 1],
                z=pos[:, 2],
                u=norm[:, 0] * -1,
                v=norm[:, 1] * -1,
                w=norm[:, 2] * -1,
                sizeref=1,
                sizemode='absolute',
                anchor='tail',
                text=labels,
                showscale=False,
                colorscale=[
                    [0, 'rgb(0, 0, 0)'],
                    [1, 'rgb(0, 0, 0)'],
                    ],
                hoverinfo='skip',
                ),
            ]
    else:
        traces = []

    traces.append(
        go.Mesh3d(
            x=pial['pos'][:, 0],
            y=pial['pos'][:, 1],
            z=pial['pos'][:, 2],
            i=pial['tri'][:, 0],
            j=pial['tri'][:, 1],
            k=pial['tri'][:, 2],
            color='pink',
            hoverinfo='skip',
            flatshading=False,
            lighting=dict(
                ambient=0.18,
                diffuse=1,
                fresnel=0.1,
                specular=1,
                roughness=0.1,
                ),
            lightposition=dict(
                x=0,
                y=0,
                z=-1,
                ),
            )
        )

    traces.append(
        go.Scatter3d(
            x=pos[:, 0],
            y=pos[:, 1],
            z=pos[:, 2],
            text=labels,
            mode='markers',
            hovertext=hovertext,
            hoverinfo='text',
            marker=marker,
            )
        )

    if angio is not None:
        traces.append(
            go.Scatter3d(
                x=angio['pos'][:, 0],
                y=angio['pos'][:, 1],
                z=angio['pos'][:, 2],
                mode='markers',
                hoverinfo='skip',
                marker=dict(
                    size=5,
                    color=log10(angio['value']),
                    symbol='diamond',
                    colorscale='Hot',
                    ),
                opacity=1,
                ))

    fig = go.Figure(
        data=traces,
        layout=go.Layout(
            showlegend=False,
            scene=dict(
                xaxis=AXIS,
                yaxis=AXIS,
                zaxis=AXIS,
                camera=dict(
                    eye=dict(
                        x=right_or_left,
                        y=0,
                        z=0.5,
                    ),
                    projection=dict(
                        type='orthographic',
                    ),
                    ),
                ),
            ),
        )

    return fig
