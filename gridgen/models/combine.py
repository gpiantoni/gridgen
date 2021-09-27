from scipy.stats import spearmanr
from numpy import array, argmax, corrcoef, arange, intersect1d, unravel_index, nansum, isnan, NaN, empty

from ..utils import normalize, match_labels, DTYPE


def merge_models(ecog, labels, vals):

    merged = empty(ecog.shape, dtype=DTYPE)
    merged.fill(NaN)
    [i0, i1] = intersect1d(merged['label'], labels, return_indices=True)[1:]
    i0r, i0c = unravel_index(i0, merged.shape)
    merged['value'][i0r, i0c] = vals[i1]

    return merged


def compare_model_with_ecog(model, ecog2d, fit):

    if model.get('functional') is None:
        chan, e, m = match_labels(ecog2d, model['morphology'])
        F = 0
        WEIGHTS = [0, ]

    else:
        chan, e, m, f = match_labels(ecog2d, model['morphology'], model['functional'])
        F = normalize(f) * fit.get('functional_weight', 1)
        if fit.get('functional_contribution') is None:
            WEIGHTS = arange(0, 110, 10)
        else:
            WEIGHTS = array(fit['functional_contribution'])

    E = normalize(e)
    M = normalize(m) * fit.get('morphology_weight', 1)

    x = []
    preds = []
    for weight in WEIGHTS:

        w = weight / 100
        prediction = w * F + (1 - w) * M
        preds.append(prediction)

        if fit['metric'] == 'parametric':
            c = corrcoef(E, prediction)[0, 1]
        else:
            c = spearmanr(E, prediction).correlation

        x.append(c)

    x = array(x)
    i = argmax(x)

    return WEIGHTS[i], x[i], chan, preds[i]


def sum_models(model, fit):

    if model.get('morphology') is None:
        func = model['functional']
        cc = nansum(func['value'])
        i = ~isnan(func['value'])

        return 100, cc, func[i]['label'], func[i]['value']

    if model.get('functional') is None:
        morph = model['morphology']
        cc = nansum(morph['value'])
        i = ~isnan(morph['value'])

        return 0, cc, morph[i]['label'], morph[i]['value']

    chan, m, f = match_labels(model['morphology'], model['functional'])
    M = m * fit['morphology_weight']
    F = f * fit['functional_weight']

    if fit.get('functional_contribution') is None:
        WEIGHTS = arange(0, 110, 10)
    else:
        WEIGHTS = array(fit['functional_contribution'])

    x = []
    preds = []
    for weight in WEIGHTS:

        w = weight / 100
        prediction = w * F + (1 - w) * M
        preds.append(prediction)
        c = nansum(prediction)
        x.append(c)

    x = array(x)
    i = argmax(x)

    return WEIGHTS[i], x[i], chan, preds[i]
