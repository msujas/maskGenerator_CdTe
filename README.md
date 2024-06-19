Scripts and programs for cosmic masking, integration and gain map correction for BM31 total scattering data on the Pilatus 2M CdTe detector. Requires a stack of images (~10 or more) to determine where the cosmics are by a statistical approach.

Run 'pip install -e .' (-e for editable, dropping will put files in site-packages, rather than reading directly) in the folder to install. This creates a GUI program called 'imagin' (Image MAsk Generator INtegrator). Run from the command prompt - choose directory, poni  file, mask file and gain map file, then run.

There are also some scripts and a Jupyter notebook which do the same. The basic script, maskGeneratorIntegraterCdTe.py, and similar Jupyter notebook will integrate a images in a single folder, the recursive script will look through all subfolders.

The azimuthalMaskGenerator.py script can be used to remove cosmics from a single image, but it's probably less accurate, and won't work well with samples with poor powder averaging. It is also very slow.

Other scrips: maskGenerator_splitting.py - for very large datasets (e.g. in situ), splitting up the mask generation into smaller subsets of images so very large arrays don't have to be stored.

![image](https://github.com/msujas/maskGenerator_CdTe/assets/79653376/9e5ec5b7-3e24-41f5-81ed-831a0ad399e1)


