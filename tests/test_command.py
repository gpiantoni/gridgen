from gridgen.bin.command import main
from gridgen.bin.parameters import REQUIRED, _JSONEncoder_path
from json import dump
from .paths import DATA_PATH, OUTPUT_PATH, ECOG_FILE

EXAMPLES = {
    "grid2d": {
        "n_rows": 3,
        "n_columns": 3,
        "direction": "TBLR",
        "chan_pattern": "chan{}"
        },
    "ecog": {
        "ecog_file": str(ECOG_FILE),
        "begtime": 0,
        "endtime": 10,
        "freq_range": [60, 90]
        },
    "mri": {
        "T1_file": str(DATA_PATH / "T1.mgz"),
        "dura_file": str(DATA_PATH / "lh_smooth.pial"),
        },
    "grid3d": {
        },
    "initial": {
        "label": "chan4",
        "RAS": [10, 42, 40],
        "rotation": 90,
        },
    }



def test_cmd():

    for cmd in list(REQUIRED)[:3]:

        params = {}
        params['output_dir'] = OUTPUT_PATH
        for grp in REQUIRED[cmd]:
            params[grp] = EXAMPLES[grp]

        if params.get('mri', {}).get('pial_file', None) is not None:
            params['morphology'] = {}

        param_json = OUTPUT_PATH / (cmd + '.json')

        with param_json.open('w') as f:
            dump(params, f, indent=2, cls=_JSONEncoder_path)

        main([str(param_json), cmd])
