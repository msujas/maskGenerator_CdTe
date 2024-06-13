import fabio
import numpy as np
#import matplotlib.pyplot as plt
import pyFAI
import os
import maskGeneratorBM31 as ifuns

file = r'C:\Users\kenneth1a\Documents\beamlineData\May2023\measurements\Si\xrd\average/average_gainCorrected.cbf'
bkgfile = r'C:\Users\kenneth1a\Documents\beamlineData\May2023\measurements\blankCapillary\0p3mm\xrd\average/average_gainCorrected.cbf'
poni = r'C:\Users\kenneth1a\Documents\beamlineData\May2023\Si00_tilt/pos05.poni'
poni = pyFAI.load(poni)
array = fabio.open(file).data
bkgarray = fabio.open(bkgfile).data

array_bkgsub = array-bkgarray
fimage = fabio.cbfimage.CbfImage(array_bkgsub)
fimage.save(file.replace('.cbf','_bkgsub.cbf'))
mask = np.where(array < 0, 1, 0)
mask = np.where(bkgarray < 0, 1, mask)
direc = os.path.dirname(file)
filename1d = f'{direc}/tmp.xye'
filename2d = f'{direc}/tmp.edf'
result = poni.integrate2d(array_bkgsub, npt_rad = 5000, npt_azim=360, mask = mask, polarization_factor=0.99, method = 'bbox', unit = '2th_deg',
                          error_model='poisson', correctSolidAngle=True, safe = False)
poni.integrate1d(array_bkgsub, filename=filename1d, npt = 5000, mask = mask, polarization_factor=0.99, method = 'bbox', unit = '2th_deg',
                          error_model='poisson', correctSolidAngle=True, safe = False)
ifuns.bubbleHeader(filename2d,*result[:3])
ifuns.clearPyFAI_header(filename1d)
#plt.figure(dpi = 300)
#plt.imshow(result[0],aspect = 'auto')
#plt.show()