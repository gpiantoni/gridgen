from gridloc.bin.command import main
from json import dump
from .paths import DATA_PATH, OUTPUT_PATH, PARAMETERS_FILE


PARAMETERS = {
    "output": str(OUTPUT_PATH),
    "grid": {
        "n_rows": 4,
        "n_columns": 5,
        "chan_pattern": "chan{}"
        },
    "ecog": {
        "file": str(DATA_PATH / "ecog.eeg"),
        "ref_chan": [],
        "begtime": 0,
        "endtime": 10,
        "freq_range": [60, 90]
        },
    "fitting": {
        "T1_file": str(DATA_PATH / "T1.mgz"),
        "dura_file": str(DATA_PATH / "lh_smooth.pial"),
        "pial_file": str(DATA_PATH / "lh.pial"),
        "angio_file": str(DATA_PATH / "angiogram.nii.gz"),
        "angio_threshold": 100,
        "initial": []
        }
    }


def test_cmd_grid2d():
    with PARAMETERS_FILE.open('w') as f:
        dump(PARAMETERS, f, indent=2)

    main([str(PARAMETERS_FILE), 'grid2d'])
