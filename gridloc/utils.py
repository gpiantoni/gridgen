from os import nice
from numpy import nanmax, nanmin, intersect1d, isnan, array


def be_nice():
    nice(10)


def normalize(x):
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
