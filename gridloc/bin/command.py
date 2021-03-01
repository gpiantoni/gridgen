"""Specify how to call gridloc from the command line"""
from pathlib import Path
from logging import getLogger, StreamHandler, Formatter, INFO, DEBUG
from argparse import ArgumentParser, RawTextHelpFormatter
from json import load, dumps
from textwrap import dedent
from numpy import set_printoptions

from plotly.offline import plot

from ..ecog.plot_ecog import plot_2d
from ..fitting import fitting
from ..matlab.pipeline import compare_to_matlab
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
    """Create the input information for the command line"""
    parser = ArgumentParser(
        description='Tools to calculate the position of ECoG grid on brain based on the neuronal activity',
        formatter_class=RawTextHelpFormatter)
    parser.add_argument('-l', '--log', default='info',
                        help='Logging level: info (default), debug')
    parser.add_argument('-o', '--output_dir',
                        help=dedent("""\
                        Output directory. Default is the current directory.
                        You can also specify it in the parameters.

                        Parameters:
                          output_dir :
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
            ecog_file : path to ECoG file
            begtime (optional) : start time in seconds from beginning of the file
            endtime (optional): end time in seconds from beginning of the file
            freq_range : low and high threshold of the frequency range of interest
            bad_channels (optional) : list of str, name of the channels to exclude

        Output:
          grid2d_ecog.tsv : values of the power spectrum per electrode
          grid2d_ecog.html : plot of the estimated activity of the power spectrum

        """))
    ecog_arg.set_defaults(function='ecog')

    # fit
    fit = list_functions.add_parser(
        'fit', help=dedent("""\
        Generate the grid, with the correct labels.

        Parameters:
          fitting :
            T1_file : path to T1 image (in particular, the T1.mgz from freesurfer)
            pial_file : path to pial surface (in particular, the lh.pial or rh.pial from freesurfer)
            dura_file : path to dura surface (for example, the smoothed pial surface)
            angio_file : path to angiogram (in NIfTI format). Optional.
            angio_threshold : value to threshold the angio_file. Optional.
            initial : start position for search
              label : label for the reference electrode
              RAS : initial location for the reference electrode
              rotation : degree of rotation of the grid (in degrees, 0° is roughly pointing up)
            correlation : 'parametric' (Pearson) or 'nonparametric' (rank)
            method : 'simplex', 'hopping', 'brute'
            range_simplex : list of 3 floats for range (x-direction, y-direction, rotation degrees)
            range_brute : list of 3 lists (for x-direction, y-direction, rotation degrees)

        Output:
          bestfit_vert{}_rot{}.label : electrode locations in freesurfer format
          bestfit_vert{}_rot{}.html : estimated activity of the best fit based on brain surface

        """))
    fit.set_defaults(function='fit')

    mat = list_functions.add_parser(
        'matlab', help=dedent("""\
        Compute values based on a conversion of the matlab code and compare the
        values with those computed by matlab

        Parameters:
          input :
            subjectInfo_file : path to subjectInfo.mat
            gridInfo_file : path to gridInfo.mat
          comparison :
            angiomap_file : path to angiomap
            prediction_file : path to predicted electrodes

        Output:
          ???
        """))
    mat.set_defaults(function='matlab')

    return parser


def main(arguments=None):
    """Main function which is called from the command line"""
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

    p_json = Path(args.parameters).resolve()
    with p_json.open() as f:
        parameters = load(f)

    if args.output_dir is not None:
        output = args.output_dir
    elif 'output_dir' in parameters:
        output = parameters['output_dir']
    else:
        output = '.'

    lg.debug(dumps(parameters, indent=2))

    parameters['output_dir'] = output
    parameters = convert_to_path(parameters)
    parameters['output_dir'].mkdir(exist_ok=True, parents=True)

    # outputs
    grid2d_tsv = parameters['output_dir'] / 'grid2d_labels.tsv'
    ecog_tsv = parameters['output_dir'] / 'grid2d_ecog.tsv'
    ecog_fig = parameters['output_dir'] / 'grid2d_ecog.html'
    transform_file = parameters['output_dir'] / 'tkras'

    if args.function == 'grid2d':
        from ..construct import make_grid_with_labels

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
            output=parameters['output_dir'],
            **parameters['fitting'])

    if args.function == 'matlab':
        compare_to_matlab(parameters)


def convert_to_path(d):
    for k, v in d.items():
        if v is None:
            continue
        elif isinstance(v, dict):
            v = convert_to_path(v)
        elif k.endswith('_file') or k.endswith('_dir'):
            d[k] = Path(v).resolve()
    return d
