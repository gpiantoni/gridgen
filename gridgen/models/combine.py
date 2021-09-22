from scipy.stats import spearmanr
from numpy import array, argmin, corrcoef, arange, intersect1d, unravel_index

from ..utils import normalize, match_labels


def merge_models(ecog, labels, vals):

    merged = ecog.copy()
    [i0, i1] = intersect1d(merged['label'], labels, return_indices=True)[1:]
    i0r, i0c = unravel_index(i0, merged.shape)
    merged['value'][i0r, i0c] = vals[i1]

    return merged


def compare_model_with_ecog(model, ecog2d, fit):

    if model['functional'] is None:
        chan, e, m = match_labels(ecog2d, model['morphology'])
        F = 0
        WEIGHTS = [0, ]

    else:
        chan, e, m, f = match_labels(ecog2d, model['morphology'], model['functional'])
        F = normalize(f)

        if fit['functional_contribution'] is None:
            WEIGHTS = arange(0, 110, 10)
        else:
            WEIGHTS = array(fit['functional_contribution'])

    E = normalize(e)
    M = normalize(m)

    if fit['morphology_weight'] == 'negative':
        E = 1 - E
    if fit['functional_weight'] == 'negative':
        F = 1 - F

    x = []
    preds = []
    for weight in WEIGHTS:

        w = weight / 100
        prediction = w * F + (1 - w) * M
        preds.append(prediction)

        if fit['correlation'] == 'parametric':
            c = corrcoef(E, prediction)[0, 1]
        else:
            c = spearmanr(E, prediction).correlation

        x.append(c)

    x = array(x)
    i = argmin(x)

    return WEIGHTS[i], x[i], chan, preds[i]
