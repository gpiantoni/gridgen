from pathlib import Path
from logging import getLogger, StreamHandler, Formatter, INFO, DEBUG
from argparse import ArgumentParser, RawTextHelpFormatter
from json import load
from textwrap import dedent
from numpy import set_printoptions

from plotly.offline import plot

from ..ecog.plot_ecog import plot_2d
from ..fitting import fitting
from ..io import (
    read_grid2d,
    write_grid2d,
    read_ecog2d,
    write_ecog2d,
    read_surface_ras_shift,
    export_transform,
    )

PKG_PATH = Path(__file__).parent

lg = getLogger('gridloc')

set_printoptions(suppress=True, precision=3)


def create_arguments():
    parser = ArgumentParser(
        description='Tools to calculate the position of ECoG grid on brain based on the neuronal activity',
        formatter_class=RawTextHelpFormatter)
    parser.add_argument('-l', '--log', default='info',
                        help='Logging level: info (default), debug')
    parser.add_argument('-o', '--output',
                        help=dedent("""\
                        Output directory. Default is the current directory.
                        You can also specify it in the parameters.

                        Parameters:
                          output :
                             path to output directory
                        """))
    parser.add_argument(
        'parameters', help=dedent("""\
        Path to file with the parameters for the analysis. The file with parameters
        should be formatted as a json file."""))
    list_functions = parser.add_subparsers(
        title='Functions',
        help='Use "all" to run all the steps')

    # grid
    grid_arg = list_functions.add_parser(
        'grid2d', help=dedent("""\
        Generate the grid, with the correct labels.

        Parameters:
          grid :
            n_rows : number of rows
            n_columns : number of columns
            direction : 'TBLR' (default), 'TBRL', 'BTLR', 'BTRL', 'LRTB', 'LRBT', 'RLTB', 'RLBT'
            chan_pattern : pattern to name the channels (it should match the naming pattern of the data)

        Output:
          grid2d_labels.tsv

        """))
    grid_arg.set_defaults(function='grid2d')

    # ecog
    ecog_arg = list_functions.add_parser(
        'ecog', help=dedent("""\
        Compute values for each electrodes based on ECoG.

        Parameters:
          ecog :
            file : path to ECoG file
            ref_chan (optional) : list of str, name of the channels to reference to
            begtime (optional) : start time in seconds from beginning of the file
            endtime (optional): end time in seconds from beginning of the file
            freq_range : low and high threshold of the frequency range of interest

        """))
    ecog_arg.set_defaults(function='ecog')

    # fit
    fit = list_functions.add_parser(
        'fit', help=dedent("""\
        Generate the grid, with the correct labels.

        Parameters:
          fitting :
            T1_file : path to T1 image (in particular, the T1.mgz from freesurfer)
            iniitial : start position for search
            pial_file : path to pial surface (in particular, the lh.pial or rh.pial from freesurfer)
            dura_file : path to dura surface (for example, the smoothed pial surface)
            angio_file : path to angiogram (in NIfTI format). Optional.
            angio_threshold : value to threshold the angio_file. Optional.

        """))
    fit.set_defaults(function='fit')

    return parser


def main(arguments=None):

    parser = create_arguments()
    args = parser.parse_args(arguments)

    # log can be info or debug
    DATE_FORMAT = '%H:%M:%S'
    if args.log[:1].lower() == 'i':
        lg.setLevel(INFO)
        FORMAT = '{asctime:<10}{message}'

    elif args.log[:1].lower() == 'd':
        lg.setLevel(DEBUG)
        FORMAT = '{asctime:<10}{levelname:<10}{filename:<40}(l. {lineno: 6d}): {message}'

    formatter = Formatter(fmt=FORMAT, datefmt=DATE_FORMAT, style='{')
    handler = StreamHandler()
    handler.setFormatter(formatter)

    lg.handlers = []
    lg.addHandler(handler)
    print(args)

    p_json = Path(args.parameters).resolve()
    with p_json.open() as f:
        parameters = load(f)

    if args.output is not None:
        output = args.output
    elif 'output' in parameters:
        output = parameters['output']
    else:
        output = '.'

    parameters['output'] = Path(output).resolve()
    parameters['output'].mkdir(exist_ok=True, parents=True)
    print(parameters)

    # outputs
    grid2d_tsv = parameters['output'] / 'grid2d_labels.tsv'
    ecog_tsv = parameters['output'] / 'grid2d_ecog.tsv'
    ecog_fig = parameters['output'] / 'grid2d_ecog.html'
    transform_file = parameters['output'] / 'tkras'

    if args.function == 'grid2d':
        from ..construct import make_grid_with_labels

        # TODO there should be an intermediate function, with more options
        grid2d = make_grid_with_labels(**parameters['grid'])
        lg.info(f'Writing labels to {grid2d_tsv}')
        write_grid2d(grid2d_tsv, grid2d)

    if args.function == 'ecog':
        from ..ecog.read_ecog import read_ecog, put_ecog_on_grid2d

        lg.info(f'Reading 2d grid from {grid2d_tsv}')
        grid2d = read_grid2d(grid2d_tsv)

        timefreq = read_ecog(**parameters['ecog'])
        ecog2d = put_ecog_on_grid2d(timefreq, grid2d)

        lg.info(f'Writing ECoG values to {ecog_tsv}')
        write_ecog2d(ecog_tsv, ecog2d)

        lg.info(f'Writing ECoG image to {ecog_fig}')
        fig = plot_2d(ecog2d, 'ecog')
        plot(fig, filename=str(ecog_fig), auto_open=False, include_plotlyjs='cdn')

    if args.function == 'fit':

        offset = read_surface_ras_shift(parameters['fitting']['T1_file'])
        export_transform(offset, transform_file)

        ecog2d = read_ecog2d(ecog_tsv, grid2d_tsv)

        fitting(
            ecog=ecog2d,
            output=parameters['output'],
            **parameters['fitting'])
