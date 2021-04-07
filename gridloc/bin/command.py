"""Specify how to call gridloc from the command line"""
from pathlib import Path
from logging import getLogger, StreamHandler, Formatter, INFO, DEBUG
from argparse import ArgumentParser, RawTextHelpFormatter
from json import load, dump
from textwrap import dedent
from numpy import set_printoptions
from datetime import datetime

from .parameters import (
    convert_to_path,
    help_template,
    JSONEncoder_path,
    prepare_template,
    validate_template,
    TEMPLATE,
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

PKG_PATH = Path(__file__).parent

lg = getLogger('gridloc')

set_printoptions(suppress=True, precision=3)


def create_arguments():
    """Create the input information for the command line"""
    parser = ArgumentParser(
        description='Tools to calculate the position of ECoG grid on brain based on the neuronal activity',
        formatter_class=RawTextHelpFormatter)
    list_functions = parser.add_subparsers(
        title='Functions')

    # grid
    subparam = list_functions.add_parser(
        'parameters', help=dedent("""\
        Generate an empty parameters json file. Fields with `null` are optional.

        Output:
          parameters.json

        """),
        description='\n'.join(help_template(TEMPLATE)),
        formatter_class=RawTextHelpFormatter,
        )
    subparam.set_defaults(function='parameters')
    subparam.add_argument(
        'parameters', help=dedent("""\
        Path to file to write with the parameters for the analysis."""))

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

    # fit
    subfun2 = list_functions.add_parser(
        'fit', help=dedent("""\
        Generate the grid, with the correct labels.

        Output:
          bestfit_vert{}_rot{}.label : electrode locations in freesurfer format
          bestfit_vert{}_rot{}.html : estimated activity of the best fit based on brain surface

        """))
    subfun2.set_defaults(function='fit')

    subfun3 = list_functions.add_parser(
        'matlab', help=dedent("""\
        Compute values based on a conversion of the matlab code and compare the
        values with those computed by matlab

        Output:
          ???
        """))
    subfun3.set_defaults(function='matlab')

    for subfun in (subfun0, subfun1, subfun2, subfun3):
        subfun.add_argument(
            '-l', '--log', default='info',
            help='Logging level: info (default), debug')
        subfun.add_argument(
            '-o', '--output_dir',
            help=dedent("""\
            Output directory. Default is the current directory.
            You can also specify it in the parameters.

            Parameters:
              output_dir :
            """))
        subfun.add_argument(
            'parameters', help=dedent("""\
            Path to file with the parameters for the analysis. The file with parameters
            should be formatted as a json file."""))

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

    p_json = Path(args.parameters).resolve()
    with p_json.open() as f:
        parameters = load(f)
    parameters = validate_template(TEMPLATE, parameters)

    if args.output_dir is not None:
        output = args.output_dir
    elif 'output_dir' in parameters:
        output = parameters['output_dir']
    else:
        output = 'gridloc_output'

    parameters['output_dir'] = output
    parameters = convert_to_path(parameters, p_json.parent)
    parameters['output_dir'].mkdir(exist_ok=True, parents=True)

    # outputs
    grid2d_tsv = parameters['output_dir'] / 'grid2d_labels.tsv'
    ecog_tsv = parameters['output_dir'] / 'grid2d_ecog.tsv'
    ecog_fig = parameters['output_dir'] / 'grid2d_ecog.html'
    transform_file = parameters['output_dir'] / 'tkras'

    if args.function in ('grid2d', ):
        from ..construct import make_grid_with_labels

        grid2d = make_grid_with_labels(**parameters['grid2d'])
        lg.info(f'Writing labels to {grid2d_tsv}')
        write_grid2d(grid2d_tsv, grid2d)

    if args.function in ('ecog', ):
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

    if args.function in ('fit', ):

        offset = read_surface_ras_shift(parameters['fit']['T1_file'])
        export_transform(offset, transform_file)

        ecog2d = read_ecog2d(ecog_tsv, grid2d_tsv)

        fitting(
            ecog=ecog2d,
            output=parameters['output_dir'],
            **parameters['fit'])

    if args.function == 'matlab':
        compare_to_matlab(parameters)

    parameters['timestamp'] = datetime.now().isoformat()
    parameters_json = parameters['output_dir'] / 'parameters.json'
    with parameters_json.open('w') as f:
        dump(parameters, f, indent=2, cls=JSONEncoder_path)
