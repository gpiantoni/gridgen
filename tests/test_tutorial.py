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


def test_ecog():
    param_json = TUTORIAL_PATH / 'ecog.json'
    main([
        str(param_json),
        '--output_dir',
        str(output),
        'ecog'
        ])


def test_grid3d():
    for param in range(7):
        param_json = TUTORIAL_PATH / f'grid3d_{param + 1}.json'
        main([
            str(param_json),
            '--output_dir',
            str(output),
            'grid3d'
            ])


def test_fit():
    param_json = TUTORIAL_PATH / 'fit.json'
    main([
        str(param_json),
        '--output_dir',
        str(output),
        'fit'
        ])

def test_fit_sum():
    param_json = TUTORIAL_PATH / 'fit_sum.json'
    main([
        str(param_json),
        '--output_dir',
        str(output),
        'fit'
        ])
