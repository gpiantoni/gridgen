## Parameters

Parameters should be saved as a `json` file.
You can pass the parameters with the command `gridloc parameters.json func` where `func` is one of the possible functions (`grid2d`, `ecog`, `fit` etc).
You can also use the command `gridloc parameters.json parameters` to generate a template `parameters.json` file.

- **output_dir**: (str, optional) path to output directory
- **grid2d** (required by command `grid2d`)
  - **n_rows**: (int) number of rows
  - **n_columns**: (int) number of columns
  - **direction**: (str) direction of the grid (wires are by definition at the bottom). Values: [TBLR, TBRL, BTLR, BTRL, LRTB, LRBT, RLTB, RLBT]. Default: TBLR
  - **chan_pattern**: (str, optional) pattern to name the channels (it should match the naming pattern of the data). Examples are "chan{}" or "chan{03d}
- **ecog** (required by command `ecog`)
  - **ecog_file**: (str) path to ECoG file
  - **freq_range**: (list) low and high threshold of the frequency range of interest. Default: [65, 95]
  - **begtime**: (float, optional) start time in seconds from the beginning of the file
  - **endtime**: (float, optional) end time in seconds from the beginning of the file
  - **bad_channels**: (list, optional) list of str, name of the channels to exclude
- **grid3d** (required by command `init`, `fit`)
  - **interelec_distance**: (float, optional) distance between the electrode centers (pitch), in mm. Default: 3
  - **maximum_angle**: (float, optional) maximum angle, in degrees, that the grid can flex, between two neighboring electrodes (elasticity of the grid). Default: 5
  - **step_angle**: (float, optional) step size, in degrees, when computing range between -`maximum_angle` and +`maximum_angle`. The smaller the step size, the faster the generation of the grid. Default: 0.25
- **mri** (required by command `init`, `fit`, `matlab`)
  - **T1_file**: (str) path to T1 image (in particular, the T1.mgz from freesurfer)
  - **pial_file**: (str) path to pial surface (in particular, the lh.pial or rh.pial from freesurfer)
  - **dura_file**: (str) path to dura surface (for example, the smoothed pial surface)
  - **angio_file**: (str, optional) path to angiogram (in NIfTI format)
  - **angio_threshold**: (float, optional) value to threshold the angio_file
- **initial** (required by command `init`, `fit`)
  - **label**: (str) label for the reference electrode
  - **RAS**: (list) initial location for the reference electrode (coordinates in MRI space)
  - **rotation**: (float) degree of rotation of the grid (in degrees, 0° is roughly pointing up)
- **fit** (required by command `fit`)
  - **method**: (str) method to use (brute includes simplex as a second step). Values: [brute, simplex]. Default: brute
  - **correlation**: (str, optional) 'parametric' (Pearson, default) or 'nonparametric' (rank). Values: [parametric, nonparametric]. Default: parametric
  - **distance**: (str, optional) how to compute the distance of the morphology. Values: [ray, minimum, view, cylinder, pdf]. Default: ray
  - **maximum_distance**: (float, optional) maximum distance between electrode and pial surface. Exact interpretation depends on the type of `morphology`. Default: 10
  - **penalty**: (float, optional) exponent when computing the penalty from the distance. Morphology = 1 / distance<sup>penalty</sup>. More simply, 1 = activity decreases linearly with distance; 2 = activity decreases with the square of the distance. Default: 1
  - **steps**
    - **x**: (float, optional) Step size in mm for x-direction, for method simplex
    - **y**: (float, optional) Step size in mm for y-direction, for method simplex
    - **rotation**: (float, optional) Step size in degrees for rotation, for method simplex
  - **ranges**
    - **x**: (list, optional) Range in mm for x-direction, in format [low, step, high], for method brute
    - **y**: (list, optional) Range in mm for y-direction, in format [low, step, high], for method brute
    - **rotation**: (list, optional) Range in degrees for rotation, in format [low, step, high], for method brute
- **matlab** (required by command `matlab`)
  - **input**
    - **subjectInfo_file**: (str) path to subjectInfo.mat (can also be subj_info.json if subjectInfo.mat is not available)
    - **gridInfo_file**: (str) path to gridInfo.mat
  - **comparison**
    - **angiomap_file**: (str, optional) path to _angiomap.mat
    - **model_file**: (str, optional) path to full models ("ROI")
    - **prediction_file**: (str, optional) path to coordinates with the best fit