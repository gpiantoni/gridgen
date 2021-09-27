from pathlib import Path
from gridgen import fitting
from .utils import simulate_data, simulate_tmap


TEST_PATH = Path(__file__).resolve().parent
ANALYSIS_PATH = TEST_PATH / 'analysis'
DATA_PATH = ANALYSIS_PATH / 'data'
TUTORIAL_PATH = TEST_PATH / 'tutorial'

SMOOTH_FILE = DATA_PATH / 'lh_smooth.pial'
PIAL_FILE = DATA_PATH / 'lh.pial'
T1_FILE = DATA_PATH / 'brain.mgz'
FUNC_FILE = DATA_PATH / 'angiogram.nii.gz'

GENERATED_PATH = ANALYSIS_PATH / 'generated'
GENERATED_PATH.mkdir(parents=True, exist_ok=True)
ECOG_FILE = GENERATED_PATH / 'ecog.eeg'
TMAP_FILE = GENERATED_PATH / 'tmap.nii.gz'
OUTPUT_PATH = ANALYSIS_PATH / 'output'
OUTPUT_PATH.mkdir(parents=True, exist_ok=True)

data = simulate_data()
data.export(ECOG_FILE, export_format='BrainVision')

simulate_tmap(T1_FILE, TMAP_FILE)

fitting.MAXITER = 5


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
