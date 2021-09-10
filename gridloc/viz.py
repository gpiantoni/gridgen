from logging import getLogger
from numpy import mean, sign, nanmin, nanmax, log10
from plotly.offline import plot
from textwrap import dedent
import plotly.graph_objects as go

from .io import WIRE
from .utils import normalize


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

    grid_file = output / 'bestfit'
    fig = plot_electrodes(pial, model['grid'], model['ecog']['ecog'], 'ecog', angio=angio)
    to_html([to_div(fig), ], grid_file)
    lg.debug(f'Exported merged model to {grid_file}')

    grid_file = output / 'morphology'
    fig0 = plot_2d(model['morpho'], 'morphology')
    fig1 = plot_electrodes(pial, model['grid'], model['morpho']['morphology'], 'morphology')
    to_html([to_div(fig0), to_div(fig1)], grid_file)

    if model['vasc'] is not None:
        grid_file = output / 'vascular'
        fig0 = plot_2d(model['vasc'], 'vasculature')
        fig1 = plot_electrodes(pial, model['grid'], model['vasc']['vasculature'], 'vasculature', angio=angio)
        to_html([to_div(fig0), to_div(fig1)], grid_file)
        lg.debug(f'Exported vascular to {grid_file}')

        merged = (model['percent_vasc'] * normalize(model['vasc']['vasculature']) + (100 - model['percent_vasc']) * normalize(model['morpho']['morphology'])) / 100
        grid_file = output / 'merged'
        fig = plot_electrodes(pial, model['grid'], merged, 'merged', angio=angio)
        to_html([to_div(fig), ], grid_file)
        lg.debug(f'Exported merged model to {grid_file}')


def plot_electrodes(pial, grid, values=None, value=None, ref_label=None, angio=None):
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

    else:

        colorbar, reversescale = default_colorbar(value)
        values = values.reshape(-1)
        marker = dict(
            size=MARKER_SIZE,
            color=values,
            colorscale='Hot',
            showscale=True,
            reversescale=reversescale,
            cmin=nanmin(values),
            cmax=nanmax(values),
            colorbar=dict(
                title=dict(
                    text=colorbar,
                    )
                ),
            )

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


def plot_2d(grid2d, value='ecog'):
    """Plot the 2D grid to a plotly html

    Parameters
    ----------
    grid2d :

    value : str
        ecog
    """
    colorbar, reversescale = default_colorbar(value)

    n_rows, n_cols = grid2d.shape
    traces = [
        go.Heatmap(
            z=grid2d[value],
            text=grid2d['label'],
            hoverinfo='text+z',
            colorscale='Hot',
            reversescale=reversescale,
            colorbar=dict(
                title=dict(
                    text=colorbar,
                    )
                ),
            ),
        ]

    annots = [{
        'x': n_cols / 2,
        'y': n_rows,
        'showarrow': False,
        'text': 'WIRES',
        'xref': 'x',
        'yref': 'y',
        'align': 'right',
        'bgcolor': 'white'
        }]

    layout = go.Layout(
        width=n_cols * 60,
        height=n_rows * 60,
        annotations=annots,
        autosize=False,
        xaxis=dict(
            visible=False,
            ),
        yaxis=dict(
            autorange='reversed',
            visible=False,
            ),
        )

    fig = go.Figure(traces, layout=layout)

    return fig


def default_colorbar(value):
    if value == 'ecog':
        reversescale = False
        colorbar = 'PSD<br>(Hz<sup>-1</sup>)'
    elif value == 'morphology':
        colorbar = 'Distance (mm)'
        reversescale = True
    elif value == 'vasculature':
        colorbar = 'Vascular<br>suppression<br>(weighted<br>angiogram<br>voxels)'
        reversescale = True
    elif value == 'merged':
        colorbar = 'estimated<br>ecog<br>activity<br>(a.u.)'
        reversescale = True

    return colorbar, reversescale


def to_div(fig):
    """Convert plotly FIG into an HTML div

    Parameters
    ----------
    fig : instance of plotly.Figure
        figure to convert

    Returns
    -------
    str
        html div, containing the figure as dynamic javascript plot
    """
    return plot(fig, output_type='div', show_link=False, include_plotlyjs=False)


def to_html(divs, filename):
    """Convert DIVs, obtained from 'to_div', into one HTML file

    Parameters
    ----------
    divs : list of divs
        list of the output of 'to_div'
    filename : path
        path of the file to write (extension should be .html). It overwrites if
        it exists
    """
    filename.parent.mkdir(exist_ok=True, parents=True)
    lg.debug(f'Saving {len(divs)} plots to {filename}')

    html = dedent('''\
        <html>
          <head>
            <meta charset="utf-8" />
          </head>
          <body>
            <div>
              <script type="text/javascript">window.PlotlyConfig = {MathJaxConfig: 'local'};</script>
              <script src="https://cdn.plot.ly/plotly-latest.min.js"></script>
        ''') + '\n'.join(divs) + dedent('''\
            </div>
          </body>
        </html>
        ''')

    with filename.with_suffix('.html').open('w') as f:
        f.write(html)
