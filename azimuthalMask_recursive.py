from maskGeneratorBM31.azimuthalMaskGenerator import runRecursive

direc = r''
ponifile = r''
maskfile = r''
gainfile = None
polarisation = 0.99
scale = 10**5
stdevs = 3
nbins = 800
ext = 'cbf'
outdir = 'xye'

runRecursive(direc, ponifile, maskfile, polarisation, gainfile, stdevs,scale, nbins, ext, outdir)
