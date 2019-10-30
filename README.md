### Create smooth surface

#### Bash
```bash
mris_fill -c -r 1 lh.pial lh_filled.mgz
mri_convert lh_filled.mgz lh_filled.nii.gz
fslmaths lh_filled -kernel 3D -dilF lh_dilated
```

#### Matlab
```matlab
addpath('/path/to/freesurfer/matlab/')
make_outer_surface(...
    'lh_dilated.nii.gz', ...
    15, ...
    'lh.outer', ...
    1)
```

#### Bash
```bash
mris_smooth -nw -n 60 lh.outer smooth.pial
```


