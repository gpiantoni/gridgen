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


def test_grid3d():

    for param in ('grid3d_4.json', ):  # ('grid3d_1.json', 'grid3d_2.json', 'grid3d_3.json', ):
        param_json = TUTORIAL_PATH / param
        main([
            str(param_json),
            '--output_dir',
            str(output),
            'grid3d'
            ])
