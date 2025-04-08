from maskGeneratorBM31.azimuthalMaskGenerator import runRecursive

direc = r''
ponifile = r''
maskfile = r''
gainfile = None
polarisation = 0.99
scale = 10**5
stdevs = 3
threshold = 100
nbins = 1500
ext = 'cbf'
outdir = 'xye_az'
cpp = False

runRecursive(direc, ponifile, maskfile, polarisation, gainfile, stdevs,scale, threshold, nbins, ext, outdir, cpp)
