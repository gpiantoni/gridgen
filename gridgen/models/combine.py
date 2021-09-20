from scipy.stats import spearmanr
from numpy import dtype, zeros, array, argmin, corrcoef, arange, NaN, intersect1d, unravel_index

from ..utils import normalize, match_labels


def merge_models(model):

    d_ = dtype([
        ('label', '<U256'),   # labels cannot be longer than 256 char
        ('value', 'f4'),
        ])
    merged = zeros((model['ecog'].shape[0], model['ecog'].shape[1]), dtype=d_)
    merged['label'] = model['ecog']['label']

    if model['functional'] is None:
        merged['value'] = normalize(model['ecog']['value'])

    else:
        merged['value'].fill(NaN)

        labels, e, m, v = match_labels(
            model['ecog'],
            model['morphology'],
            model['functional']
            )
        percent_func = model['percent_functional']
        M = normalize(m)
        V = normalize(v)
        prediction = V * percent_func / 100 + M * (100 - percent_func) / 100

        [i0, i1] = intersect1d(merged['label'], labels, return_indices=True)[1:]
        i0r, i0c = unravel_index(i0, merged.shape)
        merged['value'][i0r, i0c] = prediction[i1]

    return merged


def compare_model_with_ecog(model, ecog2d, correlation='parametric', functional_contribution=None):

    if model['functional'] is None:
        chan, e, m = match_labels(ecog2d, model['morphology'])
        F = 0
        WEIGHTS = [0, ]

    else:
        chan, e, m, f = match_labels(ecog2d, model['morphology'], model['functional'])
        F = normalize(f)

        if functional_contribution is None:
            WEIGHTS = arange(0, 110, 10)
        else:
            WEIGHTS = array(functional_contribution)

    E = normalize(e)
    M = normalize(m)

    x = []
    for weight in WEIGHTS:
        w = weight / 100
        prediction = w * F + (1 - w) * M

        if correlation == 'parametric':
            c = corrcoef(E, prediction)[0, 1]
        else:
            c = spearmanr(E, prediction).correlation

        x.append(c)

    x = array(x)
    i = argmin(x)

    return len(chan), WEIGHTS[i], x[i]
