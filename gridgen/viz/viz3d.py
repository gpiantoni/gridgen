from logging import getLogger
from numpy import mean, sign, nanmin, nanmax, log10
import plotly.graph_objects as go

from ..io import WIRE
from .utils import COLORSCALE


AXIS = dict(
    title="",
    visible=False,
    zeroline=False,
    showline=False,
    showticklabels=False,
    showgrid=False,
    )
MARKER_SIZE = 5


lg = getLogger(__name__)


def plot_electrodes(pial, grid, values=None, ref_label=None, functional=None):
    """
    """
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

        values = values['value'].reshape(-1)
        marker = dict(
            size=MARKER_SIZE,
            color=values,
            colorscale=COLORSCALE,
            showscale=True,
            cmin=nanmin(values),
            cmax=nanmax(values),
            )
        hovertext = [f'{x0}<br>{x1:0.3f}' for x0, x1 in zip(labels, values)]

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
            hovertext=hovertext,
            hoverinfo='text',
            marker=marker,
            ),
        ]

    if functional is not None:
        traces.append(
            go.Scatter3d(
                x=functional['pos'][:, 0],
                y=functional['pos'][:, 1],
                z=functional['pos'][:, 2],
                mode='markers',
                hoverinfo='skip',
                marker=dict(
                    size=5,
                    color=log10(functional['value']),
                    symbol='diamond',
                    colorscale='Hot',
                    ),
                opacity=1,
                ))

    elif False:
        """do not show Cone, it's not easy to see"""
        traces.append(
            go.Cone(
                x=pos[:, 0],
                y=pos[:, 1],
                z=pos[:, 2],
                u=norm[:, 0] * -1,
                v=norm[:, 1] * -1,
                w=norm[:, 2] * -1,
                sizeref=2,
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
            )

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
