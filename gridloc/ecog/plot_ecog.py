import plotly.graph_objs as go


def plot_ecog(ecog2d):

    n_rows, n_cols = ecog2d.shape
    traces = [
        go.Heatmap(
            z=ecog2d['ecog'],
            text=ecog2d['label'],
            hoverinfo='text+z',
            colorscale ='Hot',
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
