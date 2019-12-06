from pathlib import Path
from logging import getLogger, StreamHandler, Formatter, INFO, DEBUG
from argparse import ArgumentParser, RawTextHelpFormatter
from json import load
from textwrap import dedent

from ..io import (
    read_grid2d,
    write_grid2d,
    )
from ..construct import make_grid

PKG_PATH = Path(__file__).parent
lg = getLogger('gridloc')


def create_arguments():
    parser = ArgumentParser(
        description='Tools to compute position of ECoG grid on brain',
        formatter_class=RawTextHelpFormatter)
    parser.add_argument('-l', '--log', default='info',
                        help='Logging level: info (default), debug')
    parser.add_argument('parameters',
        help='path to file with the parameters for the analysis')
    list_functions = parser.add_subparsers(
        title='Functions',
        help='use "all" to run all the steps')

    # grid
    grid_arg = list_functions.add_parser(
        'grid2d', help=dedent("""\
        Generate the grid, with the correct labels.
        Parameters:
          grid :
            n_rows : number of rows
            n_columns : number of columns
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
            ref_chan : list of str, name of the channels to reference to
            begtime : start time in seconds from beginning of the file
            endtime : end time in seconds from beginning of the file
            chan_std_threshold : threshold of the s.d. to reject a channel
            freq_range : low and high threshold of the frequency range of interest

        """))
    ecog_arg.set_defaults(function='ecog')

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

    parameters['output'] = Path(parameters['output']).resolve()
    parameters['output'].mkdir(exist_ok=True, parents=True)
    print(parameters)

    grid2d_tsv = parameters['output'] / 'grid2d_labels.tsv'

    if args.function == 'grid2d':
        grid2d = make_grid(**parameters['grid'])
        lg.info(f'Writing labels to {grid2d_tsv}')
        write_grid2d(grid2d_tsv, grid2d)

    if args.function == 'ecog':
        pass
