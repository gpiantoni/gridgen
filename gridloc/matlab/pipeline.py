from nibabel.freesurfer.io import write_geometry
from numpy import reshape, corrcoef, isnan, mean, std
from numpy.linalg import norm
from logging import getLogger

from .vascular import calculateAngioMap
from ..io import read_ecog2d, read_surf, read_surface_ras_shift, read_elec
from .io import read_matlab
from .utils import get_initial_from_matlab
from ..fitting import fitting
from ..examine import measure_distances, measure_angles


lg = getLogger(__name__)


def compare_to_matlab(parameters):
    convert_neuralAct_surfaces(parameters)

    compare_position_in_space(parameters)

    cc = compare_ecog(parameters)
    lg.info(f'Correlation of gamma activity between matlab and python: {cc:0.3f}')

    if False:
        cc = compare_angio(parameters)
        if cc is not None:
            lg.info(f'Correlation of angiogram projection between matlab and python: {cc:0.3f}')

    compare_fitting(parameters)


def compare_fitting(parameters):

    grid2d_tsv = parameters['output_dir'] / 'grid2d_labels.tsv'
    ecog_tsv = parameters['output_dir'] / 'grid2d_ecog.tsv'
    ecog2d = read_ecog2d(ecog_tsv, grid2d_tsv)

    parameters = get_initial_from_matlab(parameters)
    parameters['fit']['method'] = 'brute'
    parameters['fit']['ranges'] = {
        'x': [-10, 2, 10],
        'y': [-10, 2, 10],
        'rotation': [-45, 5, 45],
        }
    gridv2 = fitting(
        ecog=ecog2d,
        output=parameters['output_dir'],
        **parameters['fit'])

    if parameters['matlab']['comparison']['prediction_file'] is not None:
        lg.info('Gridloc(matlab) best fit')
        gridv1 = read_elec(
            gridv2,
            parameters['matlab']['comparison']['prediction_file'])
        measure_distances(gridv1)
        measure_angles(gridv1)

        distances = norm(gridv2['pos'] - gridv1['pos'], axis=2)
        lg.info(f'Distance between best fit in v1 and best fit in v2: {mean(distances):0.3f}mm (sd {std(distances):0.3f})')


def convert_neuralAct_surfaces(parameters):

    subjectInfo = read_matlab(parameters['matlab']['input']['subjectInfo_file'])
    neuralAct = read_matlab(subjectInfo['neuralAct'])

    conversion_dir = parameters['output_dir'] / 'surfaces'
    conversion_dir.mkdir(exist_ok=True, parents=True)
    parameters['matlab']['surfaces'] = {}
    for surf_type in list(neuralAct):
        surf_file = conversion_dir / surf_type
        tri = neuralAct[surf_type]['tri'].astype(int) - 1
        write_geometry(str(surf_file), neuralAct[surf_type]['vert'], tri)
        parameters['matlab']['surfaces'][surf_type] = surf_file

    return parameters


def compare_position_in_space(parameters):
    ras_shift = read_surface_ras_shift(parameters['fit']['T1_file'])

    for surf_type in ['cortex', 'hullcortex', 'cortexcoarser']:
        surf = read_surf(parameters['matlab']['surfaces'][surf_type], ras_shift=[0, 0, 0], normals=False)
        print(f'{surf_type: >20}: [{min(surf["pos"][:, 2]): 7.3f} - {max(surf["pos"][:, 2]): 7.3f}]')

    for surf_type in ['dura_file', 'pial_file']:
        surf = read_surf(parameters['fit'][surf_type], ras_shift=ras_shift, normals=False)
        print(f'{surf_type: >20}: [{min(surf["pos"][:, 2]): 7.3f} - {max(surf["pos"][:, 2]): 7.3f}]')


def compare_ecog(parameters):

    grid2d_tsv = parameters['output_dir'] / 'grid2d_labels.tsv'
    ecog_tsv = parameters['output_dir'] / 'grid2d_ecog.tsv'
    ecog2d = read_ecog2d(ecog_tsv, grid2d_tsv)

    subjectInfo = read_matlab(parameters['matlab']['input']['subjectInfo_file'])
    gamma_mean = read_matlab(subjectInfo['gamma_mean'])

    ecog2d_orig = ecog2d.copy()
    ecog2d['ecog'] = reshape(gamma_mean, ecog2d.shape, 'F')

    i_x = ecog2d_orig['ecog'].flatten()
    i_y = ecog2d['ecog'].flatten()

    cc = corrcoef(
        i_x[~isnan(i_x) & ~isnan(i_y)],
        i_y[~isnan(i_x) & ~isnan(i_y)],
        )[0, 1]
    return cc


def compare_angio(parameters):
    subjectInfo = read_matlab(parameters['matlab']['input']['subjectInfo_file'])
    if 'angiomap_file' not in parameters['matlab']['comparison'] or parameters['matlab']['comparison']['angiomap_file'] is None:
        lg.warning('No angiomap to compare to. Skipping')
        return

    cortex = read_surf(parameters['matlab']['surfaces']['cortex'], ras_shift=[0, 0, 0], normals=True)

    [angioMap, normAngio] = calculateAngioMap(subjectInfo, subjectInfo['Tthreshold'], subjectInfo['VoxelDepth'], plotAngio=False, cortex=cortex)
    mat_angio = read_matlab(parameters['matlab']['comparison']['angiomap_file'])

    cc = corrcoef(mat_angio['angioMap'], angioMap)[0, 1]
    return cc
