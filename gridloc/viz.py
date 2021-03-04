from logging import getLogger
from numpy import mean, sign, nanmin, nanmax
from plotly.offline import plot
from textwrap import dedent
import plotly.graph_objects as go

from .io import export_grid
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


def plot_results(model, pial, ras_shift, output):

    grid_file = output / 'ecog'
    fig = plot_electrodes(pial, model['grid'], model['ecog']['ecog'])
    to_html([to_div(fig), ], grid_file)
    lg.debug(f'Exported merged model to {grid_file}')

    export_grid(model['grid'], ras_shift, grid_file)

    grid_file = output / 'morphology'
    fig0 = plot_2d(model['morpho'], 'morphology')
    fig1 = plot_electrodes(pial, model['grid'], model['morpho']['morphology'])
    to_html([to_div(fig0), to_div(fig1)], grid_file)

    if model['vasc'] is not None:
        grid_file = output / 'vascular'
        fig0 = plot_2d(model['vasc'], 'vasculature')
        fig1 = plot_electrodes(pial, model['grid'], model['vasc']['vasculature'])
        to_html([to_div(fig0), to_div(fig1)], grid_file)
        lg.debug(f'Exported vascular to {grid_file}')

        merged = (model['percent_vasc'] * normalize(model['vasc']['vasculature']) + (100 - model['percent_vasc']) * normalize(model['morpho']['morphology'])) / 100
        grid_file = output / 'merged'
        fig = plot_electrodes(pial, model['grid'], merged)
        to_html([to_div(fig), ], grid_file)
        lg.debug(f'Exported merged model to {grid_file}')


def plot_electrodes(pial, grid, values=None):
    right_or_left = sign(mean(pial['pos'][:, 0]))
    pos = grid['pos'].reshape(-1, 3)
    labels = grid['label'].reshape(-1)

    if values is None:
        marker = dict(
            size=MARKER_SIZE,
            color='black',
            )

    else:
        values = values.reshape(-1)
        marker = dict(
            size=MARKER_SIZE,
            color=values,
            colorscale='Hot',
            showscale=True,
            cmin=nanmin(values),
            cmax=nanmax(values),
            )

    traces = [
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
            ),
        go.Scatter3d(
            x=pos[:, 0],
            y=pos[:, 1],
            z=pos[:, 2],
            text=labels,
            mode='markers',
            hoverinfo='text',
            marker=marker,
            ),
        ]

    fig = go.Figure(
        data=traces,
        layout=go.Layout(
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
    """
    if not value == 'ecog':
        reversescale = True
    else:
        reversescale = False

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
                    text='PSD (Hz<sup>-1</sup>)',
                    )
                ),
            ),
        ]
    layout = go.Layout(
        width=n_cols * 60,
        height=n_rows * 60,
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
