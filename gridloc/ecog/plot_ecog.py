import plotly.graph_objs as go


def plot_2d(grid2d, value='ecog'):

    n_rows, n_cols = grid2d.shape
    traces = [
        go.Heatmap(
            z=grid2d[value],
            text=grid2d['label'],
            hoverinfo='text+z',
            colorscale='Hot',
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
