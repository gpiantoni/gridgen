"""Functions to compute relationship between electrode grid and
1) morphology (a mesh, such as pial surface)
2) functional (a volume, such as fMRI data or angiogram)
"""
from .compute import compute_model, make_grid3d_model
from .combine import merge_models, compare_model_with_ecog

__all__ = [
    'compare_model_with_ecog',
    'compute_model',
    'make_grid3d_model',
    'merge_models',
    ]
