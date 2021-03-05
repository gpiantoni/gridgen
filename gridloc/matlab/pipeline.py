"""Pipeline to compute best elec locations based on matlab code
"""
from logging import getLogger

from numpy import isnan, intersect1d, corrcoef, array

from ..io import read_surf, read_ecog2d_matlab
from .io import read_matlab
from .compat import projectElectrodes, createGrid, projectToCoarser, calculateModel, indexFuncLegacy
from ..utils import normalize


lg = getLogger(__name__)


def compute_ROI(parameters):
    lg.debug('Reading matlab and surfaces')
    subj_info = read_matlab(parameters['matlab']['input']['subjectInfo_file'])
    gridInfo = read_matlab(parameters['matlab']['input']['gridInfo_file'])
    subjstructs = {'electrodes': gridInfo['coords_ROI']['ROI_Tangent']}
    hullcortex = read_surf(parameters['matlab']['surfaces']['hullcortex'], [0, 0, 0], normals=True)

    lg.debug('Project Electrodes')
    projectedROIpoints = projectElectrodes(hullcortex, subjstructs, 25, False)

    lg.debug('Create Grid')
    ROI = createGrid(projectedROIpoints, rotation=None, turns=None, auxDims=subj_info['dims'], subj_info=subj_info, hullcortex=hullcortex)

    lg.debug('Project to Coarser')
    cortexcoarser = read_surf(parameters['matlab']['surfaces']['cortexcoarser'], [0, 0, 0], normals=True)
    ROI = projectToCoarser(ROI, cortexcoarser)

    if 'angiomap_file' not in parameters['matlab']['comparison'] or parameters['matlab']['comparison']['angiomap_file'] is None:
        normAngio = None
    else:
        angiomat = read_matlab(parameters['matlab']['comparison']['angiomap_file'])
        normAngio = angiomat['normAngio']

    lg.debug('Calculate model')
    cortex = read_surf(parameters['matlab']['surfaces']['cortex'], [0, 0, 0], normals=True)
    ROI = calculateModel(subj_info, ROI, cortex, normAngio)

    return ROI


def compute_correlations(parameters, ROI):
    lg.debug('Compute correlations')

    subj_info = read_matlab(parameters['matlab']['input']['subjectInfo_file'])
    gamma = read_ecog2d_matlab(
        subj_info['gamma_mean'],
        parameters['output_dir'] / 'grid2d_labels.tsv',
        )
    elec_label = indexFuncLegacy(subj_info)
    gamma['ecog'] = normalize(gamma['ecog'])

    ecog_label = gamma['label'].flatten(order='F')
    i_good = ~isnan(gamma['ecog'].flatten(order='F'))
    ecog_label = ecog_label[i_good]

    i_ecog, i_elec = intersect1d(ecog_label, elec_label, return_indices=True)[1:]
    gamma = gamma.flatten(order='F')[i_good][i_ecog]

    maxima = []
    for r0 in ROI:
        m = []
        for r1 in r0:
            corr = corrcoef(gamma['ecog'], r1['weights'][i_elec])[0, 1]
            r1['corr'] = corr
            m.append(corr)
        maxima.append(m)

    return array(maxima)
