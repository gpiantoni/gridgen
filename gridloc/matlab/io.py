import scipy.io as spio
from mat73 import loadmat as load73


def read_matlab(mat_file):

    try:
        out = loadold(mat_file)
        data_format = 'old'
        out = {k: v for k, v in out.items() if not k.startswith('__')}
    except NotImplementedError:
        out = load73(mat_file)
        data_format = 'h5py'

    if len(out) == 1:
        k = list(out)[0]
        val = out[k]

    else:
        return out

    if k == 'subj_info':
        for field in ('sfile', 'tfile', 'gamma_mean', 'neuralAct'):
            val[field] = (mat_file.parents[1] / val[field])

        val['dims'] = val['dims'].astype(int)

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
