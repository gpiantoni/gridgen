"""Various utilities
"""
from json import JSONEncoder
from os import nice
from numpy import nanmax, nanmin, intersect1d, isnan, array, integer
from pathlib import Path

WIRE = 'wire'


def be_nice():
    nice(10)


def normalize(x):
    """Normalize values between 0 and 1. NaNs are ignored.

    Parameters
    ----------
    x : (n,) array
        original data

    Returns
    -------
    (n,) array
        data with values between 0 and 1 (might contain NaN)
    """
    return (x - nanmin(x)) / (nanmax(x) - nanmin(x))


def match_labels(*args):
    """make sure that that the values are in the same order as the labels in
    ecog, also getting rid of bad channels"""
    args = [x.flatten('C') for x in args]

    labels = set(args[0]['label'])

    for arg in args:
        labels = labels & set(arg['label'][~isnan(arg['value'])])

    out = [list(labels), ]
    for arg in args:
        i_chan = intersect1d(arg['label'], array(list(labels)), assume_unique=False, return_indices=True)[1]
        out.append(arg['value'][i_chan])

    return out


def remove_wires(model):
    """We select rows with WIRE, so that we keep the 2d shape of the models
    """
    i_keep = model['grid']['label'] != WIRE
    i_keep = i_keep.all(axis=1)
    for field in ['ecog', 'grid', 'morphology', 'functional']:
        if model.get(field) is not None:
            model[field] = model[field][i_keep, :]

    return model


class _JSONEncoder_path(JSONEncoder):
    def default(self, obj):
        if isinstance(obj, integer):
            return int(obj)
        if isinstance(obj, Path):
            return str(obj)
