from gridgen.bin.command import main, _JSONEncoder_path
from gridgen.bin.parameters import REQUIRED
from gridgen.bin.parameters import main as main_parameters
from json import dump

from .paths import OUTPUT_PATH, EXAMPLES


def test_cmd_param():

    param_json = OUTPUT_PATH / 'template.json'
    main([str(param_json), 'parameters'])


def test_cmd():

    for cmd in ('grid2d', 'ecog'):

        params = {}
        params['output_dir'] = OUTPUT_PATH
        for grp in REQUIRED[cmd]:
            params[grp] = EXAMPLES[grp]

        param_json = OUTPUT_PATH / (cmd + '.json')

        with param_json.open('w') as f:
            dump(params, f, indent=2, cls=_JSONEncoder_path)

        main([
            '--log',
            'debug',
            str(param_json),
            cmd,
            ])


def test_help():
    main_parameters()
