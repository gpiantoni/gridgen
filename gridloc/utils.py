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

    FIELDS = 'ecog', 'morphology', 'vasculature'

    labels = set(args[0]['label'])

    for i, arg in enumerate(args):
        labels = labels & set(arg['label'][~isnan(arg[FIELDS[i]])])

    out = [list(labels), ]
    for i, arg in enumerate(args):
        i_chan = intersect1d(arg['label'], array(list(labels)), assume_unique=False, return_indices=True)[1]
        out.append(arg[FIELDS[i]][i_chan])

    return out
