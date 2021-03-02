from pathlib import Path

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
            'help': 'direction of the grid',
            },
        "chan_pattern": {
            'type': 'str',
            'necessary': False,
            'help': 'pattern to name the channels (it should match the naming pattern of the data)',
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
            },
        "begtime": {
            'type': 'float',
            'necessary': False,
            'help': 'start time in seconds from beginning of the file',
            },
        "endtime": {
            'type': 'float',
            'necessary': False,
            'help': 'end time in seconds from beginning of the file',
            },
        "bad_channels": {
            'type': 'list',
            'necessary': False,
            'help': 'list of str, name of the channels to exclude',
            },
        },
    "fit": {
        "T1_file": {
            'type': 'str',
            'necessary': True,
            'help': 'path to T1 image (in particular, the T1.mgz from freesurfer)',
            },
        "pial_file": {
            'type': 'str',
            'necessary': True,
            'help': 'path to pial surface (in particular, the lh.pial or rh.pial from freesurfer)',
            },
        "dura_file": {
            'type': 'str',
            'necessary': True,
            'help': 'path to dura surface (for example, the smoothed pial surface)',
            },
        "angio_file": {
            'type': 'str',
            'necessary': False,
            'help': 'path to angiogram (in NIfTI format)',
            },
        "angio_threshold": {
            'type': 'float',
            'necessary': False,
            'help': 'value to threshold the angio_file',
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
        "method": {
            "type": "str",
            "necessary": False,
            "values": ['simplex', 'brute', 'hopping'],
            "help": "method to use",
        },
        "correlation": {
            "type": "str",
            "necessary": False,
            "values": ['parametric', 'nonparameteric'],
            "help": "'parametric' (Pearson) or 'nonparametric' (rank)",
        },
        "ranges": {
            "x": {
                "type": "list",
                "necessary": False,
                "help": "Range in mm for x-direction. For simplex, it's a list of one float. For brute, it's [low, high] or [low, step, high]",
            },
            "y": {
                "type": "list",
                "necessary": False,
                "help": "Range in mm for y-direction. For simplex, it's a list of one float. For brute, it's [low, high] or [low, step, high]",
            },
            "rotation": {
                "type": "list",
                "necessary": False,
                "help": "Range in degrees for rotation. For simplex, it's a list of one float. For brute, it's [low, high] or [low, step, high]",
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
                out[k] = 0

            elif v['type'] == 'float':
                out[k] = 0.0

            elif v['type'] == 'str':
                out[k] = ''

            elif v['type'] == 'list':
                out[k] = []

        else:
            out[k] = prepare_template(v)
    return out


def help_template(temp):
    out = []
    for k, v in temp.items():
        if 'type' in v:
            out_ = k

            if v['type'] == 'int':
                out_ += '\t(int'

            elif v['type'] == 'float':
                out_ += '\t(float'

            elif v['type'] == 'str':
                out_ += '\t(str'

            elif v['type'] == 'list':
                out_ += '\t(list'

            if v['necessary']:
                out_ += ') '
            else:
                out_ += ', optional) '

            out_ += v['help']

            if 'values' in v:
                out_ += ' [' + ', '.join(v['values']) + ']'
            out.append(out_)

        else:
            out_ = help_template(v)
            out.append(k)
            out.extend([f'\t{x}' for x in out_])

    return out


def validate_template(temp, d):
    out = {}
    for k, v in temp.items():
        if k not in d or d[k] is None:
            continue

        if 'type' in v:
            if v['necessary'] and k not in d:
                raise ValueError(f'{k} should be defined')

            if 'values' in v and d[k] not in v['values']:
                values = ', '.join(f'{x}' for x in v['values'])
                raise ValueError(f'{k} should be one of these values: {values}')

            out[k] = d[k]

        elif k in d:
            out[k] = validate_template(v, d[k])

    return out
