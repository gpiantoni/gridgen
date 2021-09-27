# GRIDGEN

[![Python package](https://github.com/gpiantoni/gridgen/actions/workflows/python.yml/badge.svg)](https://github.com/gpiantoni/gridgen/actions/workflows/python.yml)
[![codecov](https://codecov.io/gh/gpiantoni/gridgen/branch/master/graph/badge.svg?token=6XL61XF65J)](https://codecov.io/gh/gpiantoni/gridgen)

See [documentation](https://gpiantoni.github.io/gridgen) for instructions on how to install or follow the [tutorial](https://gpiantoni.github.io/gridgen/tutorial.html).

## Quick Start

### Grid Creation
This pure-python3 package will allow you to create 3D meshes onto the convex hull (dura surface):

![grid3d](img/grid3d_1.png)

You can rotate the grid:

![grid3d rotation](img/grid3d_3.png)

or change the inter-electrode distance:

![grid3d interelec distance](img/grid3d_4.png)

### Grid Fitting
You can also use the information from the pial surface (mesh) to compute the distance to the electrodes, in different ways:

![grid3d morphology](img/grid3d_2_morpho2d.png)
![grid3d morphology](img/grid3d_2_morpho3d.png)

This information can be used to fit the ECoG values onto the most likely area of the brain where the grid was placed.

You can compute the spatial correlation between electrodes and voxels from an fMRI:

![grid3d angiogram](img/grid3d_6_scale.png)
![grid3d angiogram](img/grid3d_6.png)

and based on some fMRI activity, the script will find the location where it can record the highest amount of fMRI activity:

![fit sum](img/fit_sum.png)

## References
  - [Piantoni, G et al. *"Size of the spatial correlation between ECoG and fMRI activity."* *NeuroImage* 242(2021): 118459](https://doi.org/10.1016/j.neuroimage.2021.118459)
  - [Branco, MP et al. *"GridLoc: An automatic and unsupervised localization method for high-density ECoG grids."* *NeuroImage* 179(2018): 225-234](https://doi.org/10.1016/j.neuroimage.2018.06.050)
