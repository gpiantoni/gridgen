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
