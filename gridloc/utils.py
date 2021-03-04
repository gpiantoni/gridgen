from os import nice
from numpy import nanmax, nanmin, intersect1d


def be_nice():
    nice(10)


def normalize(x):
    return (x - nanmin(x)) / (nanmax(x) - nanmin(x))


def match_labels(ecog, *args):
    """make sure that that the values are in the same order as the labels in
    ecog, also getting rid of bad channels"""

    good = ecog['label'][ecog['good']]
    ecog_id = intersect1d(ecog['label'], good, return_indices=True)[1]
    a = ecog['ecog'].flatten('C')[ecog_id]
    output = [a, ]

    for estimate in args:
        field = (set(estimate.dtype.names) - {'label', }).pop()
        estimate_id = intersect1d(estimate['label'], good, return_indices=True)[1]
        b = estimate[field].flatten('C')[estimate_id]
        output.append(b)

    return output
