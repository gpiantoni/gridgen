from gridgen.bin.command import main, _JSONEncoder_path
from gridgen.bin.parameters import REQUIRED
from gridgen.bin.parameters import main as main_parameters
from json import dump

from .paths import OUTPUT_PATH, ECOG_FILE, T1_FILE, SMOOTH_FILE, PIAL_FILE, FUNC_FILE

EXAMPLES = {
    "grid2d": {
        "n_rows": 3,
        "n_columns": 2,
        "direction": "RLBT",
        "chan_pattern": "chan{}"
        },
    "ecog": {
        "ecog_file": ECOG_FILE,
        "begtime": 0,
        "endtime": 10,
        "freq_range": [60, 90]
        },
    "mri": {
        "T1_file": T1_FILE,
        "dura_file": SMOOTH_FILE,
        "pial_file": PIAL_FILE,
        "func_file": FUNC_FILE,
        },
    "grid3d": {
        },
    "initial": {
        "label": "chan4",
        "RAS": [-47, -1, 3],
        "rotation": 90,
        },
    "fit": {
        "method": "simplex",
        "steps": {
            "x": 0.1,
            "y": 0.1,
            "rotation": 0.1
            },
        },
    }


def test_cmd_param():

    param_json = OUTPUT_PATH / 'template.json'
    main([str(param_json), 'parameters'])


def test_cmd():

    for cmd in ('grid2d', 'ecog'):

        params = {}
        params['output_dir'] = OUTPUT_PATH
        for grp in REQUIRED[cmd]:
            params[grp] = EXAMPLES[grp]

        if params.get('mri', {}).get('pial_file', None) is not None:
            params['morphology'] = {}

        if params.get('mri', {}).get('func_file', None) is not None:
            params['functional'] = {'threshold': 100}

        param_json = OUTPUT_PATH / (cmd + '.json')

        with param_json.open('w') as f:
            dump(params, f, indent=2, cls=_JSONEncoder_path)

        main([str(param_json), cmd])


def test_help():
    main_parameters()
