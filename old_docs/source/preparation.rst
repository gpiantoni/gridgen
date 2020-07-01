Preparation of the data
=======================

For this analysis, you will need a structural MRI and an angiogram.

1. Angiogram
------------
Realign the angiogram to the structural MRI and rename it to `angiogram.nii.gz`.

Mask the angiogram with outskull, so that you don't get veins in the skin.

If the angiogram has a resolution higher than 1mm isotropic, you can also resize to 1mm.
This will speed up the calculation without impacting the results.

.. code-block:: bash

    flirt -in angiogram.nii.gz -ref angiogram.nii.gz -applyisoxfm 1 -o angio_1mm
    mv angio_1mm.nii.gz angiogram.gz

2. Freesurfer
-------------
Run `recon-all` on the structural MRI.

After running Freesurfer, you will have the files:

* `mri/brain.mgz`
* `surf/*h.surf`

3. Create a convex hull of the pial surface
-------------------------------------------

3.1 Bash
~~~~~~~~

.. code-block:: bash

   mris_fill -c -r 1 lh.pial lh_filled.mgz
   mri_convert lh_filled.mgz lh_filled.nii.gz
   fslmaths lh_filled -kernel 3D -dilF lh_dilated


3.2 Matlab
~~~~~~~~~~

.. code-block:: matlab

  addpath('/path/to/freesurfer/matlab/')
  make_outer_surface(...
      'lh_dilated.nii.gz', ...
      15, ...
      'lh.outer', ...
      1)

3.3 Bash
~~~~~~~~

.. code-block:: bash

  mris_smooth -nw -n 60 lh.outer smooth.pial
  mv lh.smooth.pial lh_smooth.pial

4. Summary of files
-------------------

* `ecog_file.ns3`
* `brain.mgz`
* `lh.pial`
* `lh_smooth.pial`
* `angiogram.nii.gz`


