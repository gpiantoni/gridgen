#!/usr/bin/env python3

from pathlib import Path
from json import load
from collections.abc import Iterable


PKG_PATH = Path(__file__).parents[2]

DIMENSIONS = 'x', 'y', 'rotation'

REQUIRED = {
    'grid2d': [
        'grid2d',
        ],
    'ecog': [
        'ecog',
        ],
    'grid3d': [
        'grid3d',
        'mri',
        'initial',
        ],
    'fit': [
        'grid3d',
        'mri',
        'initial',
        'fit',
        ],
    'matlab': [
        'mri',
        'matlab',
        ]
    }


TEMPLATE = {
    "output_dir": {
        'type': 'str',
        'necessary': False,
        'help': 'path to output directory',
        },
    "grid2d": {
        "n_rows": {
            'type': 'int',
            'necessary': True,
            'help': 'number of rows',
            },
        "n_columns": {
            'type': 'int',
            'necessary': True,
            'help': 'number of columns',
            },
        "direction": {
            'type': 'str',
            'necessary': True,
            'values': ['TBLR', 'TBRL', 'BTLR', 'BTRL', 'LRTB', 'LRBT', 'RLTB', 'RLBT'],
            'help': 'direction of the grid (wires are by definition at the bottom)',
            'default': 'TBLR'
            },
        "chan_pattern": {
            'type': 'str',
            'necessary': False,
            'help': 'pattern to name the channels (it should match the naming pattern of the data). Examples are "chan{}" or "chan{03d}',
            },
        },
    "ecog": {
        "ecog_file": {
            'type': 'str',
            'necessary': True,
            'help': 'path to ECoG file',
            },
        "freq_range": {
            'type': 'list',
            'necessary': True,
            'help': 'low and high threshold of the frequency range of interest',
            'default': [65, 95],
            },
        "begtime": {
            'type': 'float',
            'necessary': False,
            'help': 'start time in seconds from the beginning of the file',
            },
        "endtime": {
            'type': 'float',
            'necessary': False,
            'help': 'end time in seconds from the beginning of the file',
            },
        "bad_channels": {
            'type': 'list',
            'necessary': False,
            'help': 'list of str, name of the channels to exclude',
            },
        },
    "grid3d": {
        "interelec_distance": {
            "type": "float",
            "necessary": False,
            "help": "distance between the electrode centers (pitch), in mm",
            "default": 3,
            },
        "maximum_angle": {
            "type": "float",
            "necessary": False,
            "help": "maximum angle, in degrees, that the grid can flex, between two neighboring electrodes (elasticity of the grid)",
            "default": 5,
            },
        "step_angle": {
            "type": "float",
            "necessary": False,
            "help": "step size, in degrees, when computing range between -`maximum_angle` and +`maximum_angle`. The smaller the step size, the faster the generation of the grid",
            "default": 0.25,
            },
        },
    "mri": {
        "T1_file": {
            'type': 'str',
            'necessary': True,
            'help': 'path to T1 image (in particular, the T1.mgz from freesurfer). Only used to compute the mapping between T1 RAS space and surface RAS space',
            },
        "dura_file": {
            'type': 'str',
            'necessary': True,
            'help': 'path to dura surface (for example, the smoothed pial surface). This surface will be used to generate the 3D grid',
            },
        "pial_file": {
            'type': 'str',
            'necessary': False,
            'help': 'path to pial surface (in particular, the lh.pial or rh.pial from freesurfer). You need to specify key `morphology` in `parameters.json`',
            'default': None,
            },
        "func_file": {
            'type': 'str',
            'necessary': False,
            'help': 'path to angiogram or fMRI (in NIfTI format). You need to specify key `functional` in `parameters.json`',
            'default': None,
            },
        },
    "initial": {
        "label": {
            'type': 'str',
            'necessary': True,
            'help': 'label for the reference electrode',
            },
        "RAS": {
            'type': 'list',
            'necessary': True,
            'help': 'initial location for the reference electrode (coordinates in MRI space)',
            },
        "rotation": {
            'type': 'float',
            'necessary': True,
            'help': 'degree of rotation of the grid (in degrees, 0° is roughly pointing up)',
            },
        },
    "morphology": {
        "distance": {
            "type": "str",
            "necessary": False,
            "values": ['ray', 'minimum', 'view', 'cylinder', 'pdf'],
            "help": "how to compute the distance of the morphology",
            "default": "ray",
            },
        "maximum_distance": {
            'type': 'float',
            'necessary': False,
            'help': 'maximum distance between electrode and pial surface. Exact interpretation depends on the type of `morphology`',
            'default': 10,
            },
        "penalty": {
            "type": "float",
            "necessary": False,
            "help": "exponent when computing the penalty from the distance. Morphology = 1 / distance<sup>penalty</sup>. More simply, 1 = activity decreases linearly with distance; 2 = activity decreases with the square of the distance",
            "default": 1,
            },
        },
    "functional": {
        "threshold": {
            'type': 'float',
            'necessary': True,
            'help': 'value to threshold the func_file and binarize it. If None, func_file won\'t be binarized',
            'default': None,
            },
        "metric": {
            "type": "str",
            "necessary": False,
            "values": ['gaussian', 'sphere', 'inverse'],
            "help": "TODO",
            "default": "inverse",
            },
        "kernel": {
            "type": "float",
            "necessary": False,
            "help": "TODO",
            "default": 2,
            },
        },
    "fit": {
        "method": {
            "type": "str",
            "necessary": True,
            "values": ['brute', 'simplex'],
            "help": "method to use (brute includes simplex as a second step)",
            "default": "brute",
            },
        "correlation": {
            "type": "str",
            "necessary": False,
            "values": ['parametric', 'nonparametric'],
            "help": "'parametric' (Pearson, default) or 'nonparametric' (rank)",
            "default": "parametric",
            },
        "morphology_weight": {
            "type": "str",
            "necessary": False,
            "values": ['positive', 'negative'],
            "help": "TODO",
            "default": "positive",
            },
        "functional_weight": {
            "type": "str",
            "necessary": False,
            "values": ['positive', 'negative'],
            "help": "TODO",
            "default": "positive",
            },
        "functional_contribution": {
            "type": "list",
            "necessary": False,
            "help": "If both present, it's possible to combine morphology and functional models with varying weights. This parameter indicates which weights will be tested (10 = 10% of functional contribution and 90% of morphology contribution)",
            "default": [10, 20, 30, 40, 50, 60, 70, 80, 90],
            },
        "steps": {
            "x": {
                "type": "float",
                "necessary": False,
                "help": "Step size in mm for x-direction, for method simplex",
            },
            "y": {
                "type": "float",
                "necessary": False,
                "help": "Step size in mm for y-direction, for method simplex",
            },
            "rotation": {
                "type": "float",
                "necessary": False,
                "help": "Step size in degrees for rotation, for method simplex",
            },
        },
        "ranges": {
            "x": {
                "type": "list",
                "necessary": False,
                "help": "Range in mm for x-direction, in format [low, step, high], for method brute",
            },
            "y": {
                "type": "list",
                "necessary": False,
                "help": "Range in mm for y-direction, in format [low, step, high], for method brute",
            },
            "rotation": {
                "type": "list",
                "necessary": False,
                "help": "Range in degrees for rotation, in format [low, step, high], for method brute",
            },
        },
    },
    "matlab": {
        "input": {
            "subjectInfo_file": {
                'type': 'str',
                'necessary': True,
                'help': 'path to subjectInfo.mat (can also be subj_info.json if subjectInfo.mat is not available)'
                },
            "gridInfo_file": {
                'type': 'str',
                'necessary': True,
                'help': 'path to gridInfo.mat',
                },
            },
        "comparison": {
            "angiomap_file": {
                'type': 'str',
                'necessary': False,
                'help': 'path to _angiomap.mat',
                },
            "model_file": {
                'type': 'str',
                'necessary': False,
                'help': 'path to full models ("ROI")',
                },
            "prediction_file": {
                'type': 'str',
                'necessary': False,
                'help': 'path to coordinates with the best fit',
                },
            }
        }
    }


def convert_to_path(d, parent):
    """Convert strings to paths. If the path is relative, then this script makes
    it absolute in reference to "parent"

    Parameters
    ----------
    d : dict
        if a key ends in _file or _dir, it gets converted to path
    parent : path
        absolute path to use to make relative paths absolute

    Returns
    -------
    dict
        where keys ending in _file and _dir are all absolute paths
    """
    for k, v in d.items():
        if v is None:
            continue
        elif isinstance(v, dict):
            v = convert_to_path(v, parent)
        elif k.endswith('_file') or k.endswith('_dir'):
            v = Path(v)
            if not v.is_absolute():
                v = (parent / v).resolve()
            d[k] = v
    return d


def prepare_template(temp):
    out = {}
    for k, v in temp.items():
        if 'type' in v:
            if not v['necessary']:
                out[k] = None

            elif v['type'] == 'int':
                out[k] = v.get('default', 0)

            elif v['type'] == 'float':
                out[k] = v.get('default', 0.0)

            elif v['type'] == 'str':
                out[k] = v.get('default', '')

            elif v['type'] == 'list':
                out[k] = v.get('default', [])

        else:
            out[k] = prepare_template(v)
    return out


def help_template(temp):
    modules = _invert_dict(REQUIRED)

    out = []
    for k, v in temp.items():
        if 'type' in v:
            out_ = f'- **{k}**: '

            if v['type'] == 'int':
                out_ += '(int'

            elif v['type'] == 'float':
                out_ += '(float'

            elif v['type'] == 'str':
                out_ += '(str'

            elif v['type'] == 'list':
                out_ += '(list'

            if v['necessary']:
                out_ += ') '
            else:
                out_ += ', optional) '

            out_ += v['help']

            if 'values' in v:
                out_ += '. Values: [' + ', '.join(v['values']) + ']'

            if 'default' in v:
                out_ += f'. Default: {v["default"]}'

            out.append(out_)

        else:
            out_ = help_template(v)
            if k in modules:
                extra = ' (required by command ' + ', '.join(f'`{x}`' for x in modules[k]) + ')'
            else:
                extra = ''
            out.append(f'- **{k}**' + extra)
            out.extend([f'  {x}' for x in out_])

    return out


def validate_template(temp, d):
    out = {}
    for k, v in temp.items():

        if 'method' in d:
            if d['method'] == 'brute':
                for dim in DIMENSIONS:
                    if d.get('ranges', {}).get(dim, None) is None:
                        raise ValueError(f'{dim} should be defined in "ranges" with method "brute"')

                d.pop('steps', None)

            if d['method'] == 'simplex':
                for dim in DIMENSIONS:
                    if d.get('steps', {}).get(dim, None) is None:
                        raise ValueError(f'{dim} should be defined in "steps" with method "simplex"')

                d.pop('ranges', None)

        if k not in d or d[k] is None:
            if v.get('default', False):
                out[k] = v['default']

        elif 'type' in v:
            if v['necessary']:
                if k not in d:
                    raise ValueError(f'{k} should be defined')
                if isinstance(d[k], Iterable) and len(d[k]) == 0:
                    raise ValueError(f'{k} should be defined')

            if 'values' in v and d[k] not in v['values']:
                values = ', '.join(f'{x}' for x in v['values'])
                raise ValueError(f'{k} should be one of these values: {values}')

            out[k] = d[k]

        elif k in d:
            out[k] = validate_template(v, d[k])

    return out


def _invert_dict(d):
    new_dict = {}
    for k, v in d.items():
        for x in v:
            new_dict.setdefault(x, []).append(k)
    return new_dict


def parse_parameters(parameters, function, output_dir=None):
    p_json = Path(parameters).resolve()
    with p_json.open() as f:
        parameters = load(f)

    required = REQUIRED[function].copy()
    if 'mri' in parameters:
        if parameters['mri'].get('func_file', None) is not None:
            required.append('functional')
        if parameters['mri'].get('pial_file', None) is not None:
            required.append('morphology')

    for k in required:
        if k not in parameters:
            raise ValueError(f'You need to specify "{k}" when running {function}')
        else:
            parameters[k] = validate_template(TEMPLATE[k], parameters[k])

    if output_dir is not None:
        output = output_dir
    elif 'output_dir' in parameters:
        output = parameters['output_dir']
    else:
        output = 'gridgen_output'

    parameters['output_dir'] = output
    parameters = convert_to_path(parameters, p_json.parent)
    parameters['output_dir'].mkdir(exist_ok=True, parents=True)

    # move one parameter to mri because it's more straightforward
    if parameters.get('mri', {}).get('func_file', None) is not None:
        parameters['mri']['func_threshold'] = parameters['functional'].pop('threshold', None)

    return parameters


if __name__ == '__main__':

    md = [
        '## Parameters',
        '',
        'Parameters should be saved as a `json` file.',
        'You can pass the parameters with the command `gridgen parameters.json func` where `func` is one of the possible functions (`grid2d`, `ecog`, `fit` etc).',
        'You can also use the command `gridgen parameters.json parameters` to generate a template `parameters.json` file.',
        '',
        ]
    md = md + help_template(TEMPLATE)
    parameters_file = PKG_PATH / 'docs' / 'parameters.md'

    with parameters_file.open('w') as f:
        f.write('\n'.join(md))
