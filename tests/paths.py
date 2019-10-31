from pathlib import Path

TEST_PATH = Path(__file__).resolve().parent
DATA_PATH = TEST_PATH / 'data'
GENERATORS_PATH = DATA_PATH / 'generators'
GENERATORS_PATH.mkdir(exist_ok=True)
