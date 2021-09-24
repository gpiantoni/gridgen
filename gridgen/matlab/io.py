import scipy.io as spio
from json import load as json_load
from numpy import array, c_, isnan, zeros, genfromtxt, intersect1d, unravel_index
try:
    from mat73 import loadmat as load73
except ModuleNotFoundError:
    load73 = None

from ..utils import DTYPE, DTYPE_ECOG


def read_ecog2d_matlab(gamma_file, grid_file):
    """Read the values of ECoG analysis

    Parameters
    ----------
    ecog_file : Path
        file with ecog data from matlab ('gamma mean')
    grid_file : Path
        file with labels (in 2d)

    Returns
    -------
    ecog2d : 2d ndarray
        ecog (n_rows, n_columns) with fields (label, ecog)
    """
    from ..io import read_grid2d

    grid2d = read_grid2d(grid_file)
    gamma_mean = read_matlab(gamma_file)

    ecog_on_grid = zeros(grid2d.shape, dtype=DTYPE_ECOG)
    ecog_on_grid['label'] = grid2d['label']
    ecog_on_grid['value'] = gamma_mean.reshape(grid2d.shape, order='F')
    ecog_on_grid['good'] = ~isnan(ecog_on_grid['value'])

    return ecog_on_grid


def read_elec(grid2d, elec_file):
    """Read electrode locations and match them to a grid2d

    Parameters
    ----------
    grid2d : instance of grid2d
        grid2d with labels
    elec_file : path to .mat or .tsv

    Returns
    -------
    instance of grid2d
        where 'pos' are taken from elec_file
    """
    grid2d = grid2d.copy()  # prevents changing the input file

    if elec_file.suffix == '.mat':
        from .matlab.io import read_matlab

        xyz = read_matlab(elec_file)
        labels = array([f'chan{x + 1}' for x in range(xyz.shape[0])])

    elif elec_file.suffix == '.tsv':
        elec = genfromtxt(elec_file, skip_header=1, dtype=DTYPE)
        labels = elec['name']
        xyz = c_[elec['x'], elec['y'], elec['z']]

    i_grid, i_mat = intersect1d(grid2d['label'], labels, return_indices=True)[1:]
    i = unravel_index(i_grid, grid2d.shape)
    grid2d['pos'][i] = xyz[i_mat]

    return grid2d


def read_matlab(mat_file):
    """Read either matlab or json file with the same info.

    Notes
    -----
    1) If it's mat file with coordsPred, it returns only this variable and
    discards the other ones.

    2) it can also reas json files. This is useful when subjectInfo.mat was not
    created or has an error. So create a subj_info.json file with the fields:
    {
        "dims": [8, 16],
        "intElec": [3, 3],
        "hemiVect": {
            "hemi": "r",
            "side": "d"
            },
        "gamma_mean": "",
        "neuralAct": "",
        "sfile": "",
        "tfile": "",
        "Tthreshold": 50,
        "VoxelDepth": 8
    }
    """
    if mat_file.suffix == '.json':
        with mat_file.open() as f:
            out = {mat_file.stem: json_load(f)}

    else:

        try:
            out = loadold(mat_file)
            out = {k: v for k, v in out.items() if not k.startswith('__')}
        except NotImplementedError:
            if load73 is None:
                raise ModuleNotFoundError('You need to install "mat73"')
            out = load73(mat_file)

    if len(out) == 1:
        k = list(out)[0]
        val = out[k]

    elif 'coordsPred' in out:
        return out['coordsPred']

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
