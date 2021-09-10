from logging import getLogger
from plotly.offline import plot
from textwrap import dedent

lg = getLogger(__name__)


def default_colorbar(value):
    if value == 'ecog':
        reversescale = False
        colorbar = 'PSD (Hz<sup>-1</sup>)'
    elif value == 'morphology':
        colorbar = 'Distance (mm)'
        reversescale = True
    elif value == 'vasculature':
        colorbar = 'Vascular suppression (weighted average of neighboring voxels)'
        reversescale = True
    elif value == 'merged':
        colorbar = 'estimated ecog activity (a.u.)'
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
