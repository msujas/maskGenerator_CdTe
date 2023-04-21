A Jupyter Notebook for automatically generating masks for cosmic radation on the Pilatus 2M CdTe detector on BM31, SNBL, ESRF. Requires a series of images (~20) and filters by the standard deviation on each pixel. The main script is maskGeneratorIntegraterCdTe.ipynb. maskGeneratorCdTe_recursive.py is similar and will run through all folders from a specified directory.

The azimuthalMaskGenerator notebook can be used to remove cosmics from a single image, but it's probably less accurate, and won't work well with samples with poor powder averaging. It is also very slow.

Other scrips: maskGenerator_splitting.py - for very large datasets, splitting up the mask generation into smaller subsets of images so very large arrays don't have to be stored.
