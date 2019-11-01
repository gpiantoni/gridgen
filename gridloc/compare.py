from numpy.ma import masked_invalid, corrcoef
from numpy import arange, ptp, array

from gridloc.morphology.distance import compute_distance
from gridloc.vascular.sphere import vascular_model


def compare_models(grid, pial, angio, offset, gamma):

    m_morpho = compute_distance(grid, pial)
    m_vasc = vascular_model(grid, angio, -1 * offset)

    x = []
    for weight in arange(0.1, 1, 0.1):
        prediction = weight * normalize(m_vasc) + (1 - weight) * normalize(m_morpho)
        x.append(corrcoef_nan(-1 * prediction, gamma))

    return array(x).min()


def corrcoef_nan(A, B):
    a = masked_invalid(A)
    b = masked_invalid(B)

    msk = (~a.mask & ~b.mask)
    return corrcoef(a[msk], b[msk])[0, 1]


def normalize(x):
    return (x - x.min()) / ptp(x)
