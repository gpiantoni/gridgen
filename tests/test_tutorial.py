from .paths import TUTORIAL_PATH
from gridgen.bin.command import main


output = str(TUTORIAL_PATH / 'output')

def test_grid2d():
    param_json = TUTORIAL_PATH / 'grid2d.json'
    main([
        str(param_json),
        '--output_dir',
        str(output),
        'grid2d'
        ])
