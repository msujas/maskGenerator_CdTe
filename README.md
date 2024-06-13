Run 'pip install -e .' in the folder to install. This creates a GUI program called 'imagin' (Image MAsk Generator INtegrator). Run from the command prompt - choose directory, poni  file, mask file and gain map file, then run.

A Jupyter Notebook for automatically generating masks for cosmic radation on the Pilatus 2M CdTe detector on BM31, SNBL, ESRF. Requires a series of images (~20) and filters by the standard deviation on each pixel. The main script is maskGeneratorIntegraterCdTe.ipynb. maskGeneratorCdTe_recursive.py is similar and will run through all folders from a specified directory.

The azimuthalMaskGenerator.py script can be used to remove cosmics from a single image, but it's probably less accurate, and won't work well with samples with poor powder averaging. It is also very slow.

Other scrips: maskGenerator_splitting.py - for very large datasets, splitting up the mask generation into smaller subsets of images so very large arrays don't have to be stored.

![image](https://github.com/msujas/maskGenerator_CdTe/assets/79653376/9e5ec5b7-3e24-41f5-81ed-831a0ad399e1)


