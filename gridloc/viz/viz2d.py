from logging import getLogger
import plotly.graph_objects as go

from .utils import to_div, default_colorbar
from ..utils import match_labels, normalize


lg = getLogger(__name__)


def plot_grid2d(grid2d, value='ecog'):
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
                    text=colorbar.replace(' ', '<br>'),
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


def plot_scatter(model):

    labels, e, m, v = match_labels(
        model['ecog'],
        model['morpho'],
        model['vasc']
        )

    divs = []

    fig = plot_correlation(labels, e, m, 'ecog', 'morphology')
    divs.append(to_div(fig))

    fig = plot_correlation(labels, e, v, 'ecog', 'vasculature')
    divs.append(to_div(fig))

    fig = plot_prediction(labels, e, m, v, percent_vasc=model['percent_vasc'])
    divs.append(to_div(fig))

    return divs


def plot_correlation(labels, x, y, xname, yname):

    traces = [
        go.Scatter(
            x=x,
            y=y,
            mode='markers',
            text=labels,
            marker=dict(
                color='black',
            ),
        )
        ]

    layout = go.Layout(
        height=700,
        width=700,
        title=f'Correlation between {xname} and {yname} (not normalized)',
        xaxis=dict(
            title=default_colorbar(xname)[0],
            ),
        yaxis=dict(
            title=default_colorbar(yname)[0],
            ),
        )

    return go.Figure(data=traces, layout=layout)


def plot_prediction(labels, e, m, v=None, percent_vasc=None):

    E = normalize(e)
    M = normalize(m)
    title = 'Correlation between ecog and prediction (normalized)'
    if v is not None:
        V = normalize(v)
        prediction = V * percent_vasc / 100 + M * (100 - percent_vasc) / 100
        title = title + f'<br>{100 - percent_vasc:.0f}% morphology, {percent_vasc:.0f}% vasculature'
    else:
        prediction = M

    traces = [
        go.Scatter(
            y=1 - 1 * prediction,
            x=E,
            mode='markers',
            text=labels,
            marker=dict(
                color='black',
            ),
        )
        ]

    layout = go.Layout(
        title=title,
        height=700,
        width=700,
        yaxis=dict(
            title='predicted values',
            range=(-0.05, 1.05),
            ),
        xaxis=dict(
            title='normalized observed values',
            range=(-0.05, 1.05),
            ),
        )

    return go.Figure(data=traces, layout=layout)
