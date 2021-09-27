"""Specify how to call gridgen from the command line"""
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
    )

from ..fitting import fitting
from ..matlab import compare_to_matlab
from ..models import make_grid3d_model
from ..viz import to_html, to_div, plot_grid2d, plot_grid3d
from ..grid2d import make_grid_with_labels
from ..ecog import read_ecog, put_ecog_on_grid2d
from ..utils import _JSONEncoder_path, remove_wires
from ..io import (
    read_mri,
    read_grid2d,
    write_grid2d,
    read_ecog2d,
    write_ecog2d,
    read_surface_ras_shift,
    export_transform,
    export_electrodes,
    )

lg = getLogger('gridgen')

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

    subfun2 = list_functions.add_parser(
        'grid3d', help=dedent("""\
        Generate 3D grid.

        Input:
          grid2d_labels.tsv (from grid2d)

        Output (in subfolder):
          electrodes.tsv : electrode locations in T1 space
          electrodes.label : electrode locations for freeview
          electrodes.fcsv : electrode locations for 3DSlicer
          electrodes.html : interactive plot with electrode locations

        """))
    subfun2.set_defaults(function='grid3d')

    # ecog
    subfun1 = list_functions.add_parser(
        'ecog', help=dedent("""\
        Compute values for each electrodes based on ECoG.

        Input:
          grid2d_labels.tsv (from grid2d)

        Output:
          grid2d_ecog.tsv : values of the power spectrum per electrode
          grid2d_ecog.html : plot of the estimated activity of the power spectrum

        """))
    subfun1.set_defaults(function='ecog')

    subfun3 = list_functions.add_parser(
        'fit', help=dedent("""\
        Fit the ecog values to the surface

        Output:
          parameters.json : summary of the parameters used for the fit
          results.json : summary of the results to recompute the grid and fit
          electrodes.tsv : electrode locations for the best fit in T1 space
          electrodes.label : electrode locations for the best fit for freeview
          electrodes.fcsv : electrode locations for the best fit for 3DSlicer
          projected.html : interactive plot with electrode locations for the best fit (with ECoG values)
          morphology.html : the values for each electrode based on morphology of the pial surface
          functional.html : the values for each electrode based on functional MRI
          merged.html : the combined values of functional and morphology values

        """))
    subfun3.set_defaults(function='fit')

    subfun4 = list_functions.add_parser(
        'matlab', help=dedent("""\
        Compute values based on a conversion of the matlab code and compare the
        values with those computed by matlab
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

    start_time = datetime.now()

    if args.function in ('ecog', 'grid3d'):
        lg.info(f'Reading 2d grid from {grid2d_tsv}')
        grid2d = read_grid2d(grid2d_tsv)

    if args.function in ('grid3d', 'fit'):
        offset = read_surface_ras_shift(parameters['mri']['T1_file'])
        export_transform(offset, transform_file)
        mris = read_mri(**parameters['mri'])

    if args.function == 'grid2d':
        grid2d = make_grid_with_labels(**parameters['grid2d'])
        lg.info(f'Writing labels to {grid2d_tsv}')
        write_grid2d(grid2d_tsv, grid2d)

    if args.function == 'ecog':
        timefreq = read_ecog(**parameters['ecog'])
        ecog2d = put_ecog_on_grid2d(timefreq, grid2d)

        lg.info(f'Writing ECoG values to {ecog_tsv}')
        write_ecog2d(ecog_tsv, ecog2d)

        lg.info(f'Writing ECoG image to {ecog_fig}')
        fig = plot_grid2d(ecog2d)
        to_html([to_div(fig), ], ecog_fig)

    if args.function == 'grid3d':

        grid2d = read_grid2d(grid2d_tsv)
        output_dir = parameters['output_dir'] / ('grid3d_' + start_time.strftime('%Y%m%d_%H%M%S'))
        output_dir.mkdir(parents=True)
        lg.info(f'Writing grid3d to {output_dir}')
        parameters['output_dir'] = output_dir

        model = make_grid3d_model(
            grid2d=grid2d,
            mris=mris,
            grid3d=parameters['grid3d'],
            initial=parameters['initial'],
            morphology=parameters.get('morphology', {}),
            functional=parameters.get('functional', {}),
            )

        plot_grid3d(parameters, mris, model)
        model = remove_wires(model)
        export_electrodes(output_dir, model, mris)

    if args.function == 'fit':
        if parameters['fit']['metric'] == 'sum':
            ecog2d = read_grid2d(grid2d_tsv)
        else:
            ecog2d = read_ecog2d(ecog_tsv, grid2d_tsv)

        folder_name = '_'.join(str(parameters['fit'][k]) for k in ('metric', 'method'))
        output_dir = parameters['output_dir'] / ('fit_' + start_time.strftime('%Y%m%d_%H%M%S') + '_' + folder_name)
        output_dir.mkdir(parents=True)
        lg.info(f'Writing fitting results to {output_dir}')

        parameters['timestamp'] = start_time.isoformat()
        parameters_json = output_dir / 'parameters.json'
        with parameters_json.open('w') as f:
            dump(parameters, f, indent=2, cls=_JSONEncoder_path)

        fitting(
            output=output_dir,
            ecog=ecog2d,
            mris=mris,
            grid3d=parameters['grid3d'],
            initial=parameters['initial'],
            fit=parameters['fit'],
            morphology=parameters.get('morphology', {}),
            functional=parameters.get('functional', {}),
            )

    if args.function == 'matlab':
        compare_to_matlab(parameters)
