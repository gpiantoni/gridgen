"""Specify how to call gridloc from the command line"""
from pathlib import Path
from logging import getLogger, StreamHandler, Formatter, INFO, DEBUG
from argparse import ArgumentParser, RawTextHelpFormatter
from json import dump
from textwrap import dedent
from numpy import set_printoptions
from datetime import datetime

from .parameters import (
    prepare_template,
    parse_parameters,
    TEMPLATE,
    _JSONEncoder_path
    )
from ..fitting import fitting
from ..matlab.comparison import compare_to_matlab
from ..viz import to_html, to_div, plot_2d
from ..io import (
    read_grid2d,
    write_grid2d,
    read_ecog2d,
    write_ecog2d,
    read_surface_ras_shift,
    export_transform,
    )

lg = getLogger('gridloc')

set_printoptions(suppress=True, precision=3)


def create_arguments():
    """Create the input information for the command line"""
    parser = ArgumentParser(
        description='Tools to calculate the position of ECoG grid on brain based on the neuronal activity',
        formatter_class=RawTextHelpFormatter)

    parser.add_argument(
        'parameters', help=dedent("""\
        Path to file with the parameters for the analysis. The file with parameters
        should be formatted as a json file."""))
    list_functions = parser.add_subparsers(
        title='Functions')
    parser.add_argument(
        '-l', '--log', default='info',
        help='Logging level: info (default), debug')
    parser.add_argument(
        '-o', '--output_dir',
        help=dedent("""\
        Output directory. Default is the current directory.
        You can also specify it in the parameters.

        Parameters:
          output_dir :
        """))

    # create parameters
    subparam = list_functions.add_parser(
        'parameters', help=dedent("""\
        Generate an empty parameters json file. Fields with `null` are optional.

        Output:
          parameters.json

        """),
        )
    subparam.set_defaults(function='parameters')

    # grid
    subfun0 = list_functions.add_parser(
        'grid2d', help=dedent("""\
        Generate the grid, with the correct labels.

        Output:
          grid2d_labels.tsv

        """))
    subfun0.set_defaults(function='grid2d')

    # ecog
    subfun1 = list_functions.add_parser(
        'ecog', help=dedent("""\
        Compute values for each electrodes based on ECoG.

        Output:
          grid2d_ecog.tsv : values of the power spectrum per electrode
          grid2d_ecog.html : plot of the estimated activity of the power spectrum

        """))
    subfun1.set_defaults(function='ecog')

    subfun2 = list_functions.add_parser(
        'init', help=dedent("""\
        Fit the ecog values to the surface

        Output:
            ???

        """))
    subfun2.set_defaults(function='init')

    subfun3 = list_functions.add_parser(
        'fit', help=dedent("""\
        Fit the ecog values to the surface

        Output:
            ???

        """))
    subfun3.set_defaults(function='fit')

    subfun4 = list_functions.add_parser(
        'matlab', help=dedent("""\
        Compute values based on a conversion of the matlab code and compare the
        values with those computed by matlab

        Output:
          ???
        """))
    subfun4.set_defaults(function='matlab')

    return parser


def main(arguments=None):
    """Main function which is called from the command line"""
    parser = create_arguments()
    args = parser.parse_args(arguments)

    if args.function == 'parameters':
        parameters = prepare_template(TEMPLATE)
        p_json = Path(args.parameters).resolve().with_suffix('.json')
        with p_json.open('w') as f:
            dump(parameters, f, indent=2)
        return

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

    parameters = parse_parameters(args.parameters, args.function, args.output_dir)

    # outputs
    grid2d_tsv = parameters['output_dir'] / 'grid2d_labels.tsv'
    ecog_tsv = parameters['output_dir'] / 'grid2d_ecog.tsv'
    ecog_fig = parameters['output_dir'] / 'grid2d_ecog.html'
    transform_file = parameters['output_dir'] / 'tkras'

    if args.function == 'grid2d':
        from ..construct import make_grid_with_labels

        grid2d = make_grid_with_labels(**parameters['grid2d'])
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
        to_html([to_div(fig), ], ecog_fig)

    if args.function in ('init', 'fit'):

        offset = read_surface_ras_shift(parameters['mri']['T1_file'])
        export_transform(offset, transform_file)

        ecog2d = read_ecog2d(ecog_tsv, grid2d_tsv)

        if args.function == 'fit':
            start_time = datetime.now()
            output_dir = parameters['output_dir'] / ('bestfit_' + parameters['fit']['method'] + '_' + parameters['fit']['correlation'] + '_' + start_time.strftime('%Y%m%d_%H%M%S'))
            output_dir.mkdir(parents=True)

            parameters['timestamp'] = start_time.isoformat()
            parameters_json = output_dir / 'parameters.json'
            with parameters_json.open('w') as f:
                dump(parameters, f, indent=2, cls=_JSONEncoder_path)
        else:
            output_dir = parameters['output_dir']

        fitting(
            ecog=ecog2d,
            output=output_dir,
            init=args.function == 'init',
            **parameters['fit'])

    if args.function == 'matlab':
        compare_to_matlab(parameters)
