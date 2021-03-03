import scipy.io as spio
from mat73 import loadmat as load73
from json import load as json_load
from numpy import array


def read_matlab(mat_file):
    """Read either matlab or json file with the same info.

    Json file is useful when subjectInfo.mat was not created. So create a
    subj_info.json file with the fields:
        'dims': [8, 16],
        'intElec': [3, 3],
        'hemiVect': {
            'hemi': 'r',
            'side': 'd',
            },
        'gamma_mean': '',
        'neuralAct': '',
        'sfile': '',
        'tfile': '',
        'Tthreshold': 50,
        'VoxelDepth': 8,
    """

    if mat_file.suffix == '.json':
        with mat_file.open() as f:
            out = {mat_file.stem: json_load(f)}

    else:

        try:
            out = loadold(mat_file)
            out = {k: v for k, v in out.items() if not k.startswith('__')}
        except NotImplementedError:
            out = load73(mat_file)

    if len(out) == 1:
        k = list(out)[0]
        val = out[k]

    else:
        return out

    if k == 'subj_info':
        for field in ('sfile', 'tfile', 'gamma_mean', 'neuralAct'):
            if val[field] is None or val[field] == '':
                val[field] = None
            else:
                if val[field].startswith('.'):
                    val[field] = '.' + val[field]
                val[field] = (mat_file.parent / val[field])

        val['dims'] = array(val['dims']).astype(int)

    return val


def loadold(filename):
    '''
    https://stackoverflow.com/a/8832212

    this function should be called instead of direct spio.loadmat
    as it cures the problem of not properly recovering python dictionaries
    from mat files. It calls the function check keys to cure all entries
    which are still mat-objects
    '''
    data = spio.loadmat(filename, struct_as_record=False, squeeze_me=True)
    return _check_keys(data)

def _check_keys(dict):
    '''
    checks if entries in dictionary are mat-objects. If yes
    todict is called to change them to nested dictionaries
    '''
    for key in dict:
        if isinstance(dict[key], spio.matlab.mio5_params.mat_struct):
            dict[key] = _todict(dict[key])
    return dict

def _todict(matobj):
    '''
    A recursive function which constructs from matobjects nested dictionaries
    '''
    dict = {}
    for strg in matobj._fieldnames:
        elem = matobj.__dict__[strg]
        if isinstance(elem, spio.matlab.mio5_params.mat_struct):
            dict[strg] = _todict(elem)
        else:
            dict[strg] = elem
    return dict
