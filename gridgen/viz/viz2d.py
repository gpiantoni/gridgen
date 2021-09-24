from logging import getLogger
import plotly.graph_objects as go

from .utils import to_div, COLORSCALE
from ..utils import match_labels


lg = getLogger(__name__)

modalities = 'morphology', 'functional'


def plot_grid2d(grid2d):
    """Plot the 2D grid to a plotly html

    Parameters
    ----------
    grid2d :

    """
    n_rows, n_cols = grid2d.shape
    traces = [
        go.Heatmap(
            z=grid2d['value'],
            text=grid2d['label'],
            hoverinfo='text+z',
            colorscale=COLORSCALE,
            ),
        ]

    annots = [{
        'x': n_cols / 2 - 0.5,
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


def plot_scatter(model, fit):

    divs = []

    for mod in modalities:
        if model.get(mod, None) is None:
            continue
        fig = plot_correlation(model, 'ecog', mod)

        if fit[f'{mod}_weight'] == 'negative':
            fig.update_layout(yaxis=dict(autorange='reversed'))
        divs.append(to_div(fig))

    fig = plot_correlation(model, 'ecog', 'merged')
    divs.append(to_div(fig))

    return divs


def plot_correlation(model, xname, yname):

    labels, x, y = match_labels(
        model[xname],
        model[yname],
        )

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
        title=f'Correlation between {xname} and {yname}',
        xaxis=dict(
            title=xname,
            ),
        yaxis=dict(
            title=yname,
            ),
        )

    return go.Figure(data=traces, layout=layout)
