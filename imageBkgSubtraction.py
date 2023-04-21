import fabio
import numpy as np
#import matplotlib.pyplot as plt
import pyFAI
import os
import integrationFunctions as ifuns
direc = r'C:\Users\kenneth1a\Documents\beamlineData\a311189_pdf_2/'
os.chdir(direc)
file = r'39_glycine_B1mm\pdf\xrd\average/average_gainCorrected.cbf'
bkgfile = r'44_empty_B1mm\pdf\xrd\average/average_gainCorrected.cbf'
poni = 'Si_0_15tilt.poni'
poni = pyFAI.load(poni)
array = fabio.open(file).data
bkgarray = fabio.open(bkgfile).data

array_bkgsub = array-bkgarray
mask = np.where(array < 0, 1, 0)
mask = np.where(bkgarray < 0, 1, mask)

result = poni.integrate2d(array_bkgsub, npt_rad = 5000, npt_azim=360, mask = mask, polarization_factor=0.99, method = 'bbox', unit = '2th_deg',
                          error_model='poisson', correctSolidAngle=True)
poni.integrate1d(array_bkgsub, filename='tmp.xye', npt = 5000, mask = mask, polarization_factor=0.99, method = 'bbox', unit = '2th_deg',
                          error_model='poisson', correctSolidAngle=True)
ifuns.bubbleHeader('tmp.edf',*result[:3])
#plt.figure(dpi = 300)
#plt.imshow(result[0],aspect = 'auto')
#plt.show()