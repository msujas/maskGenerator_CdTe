import numpy as np
import fabio
import matplotlib.pyplot as plt
import pyFAI
import pyFAI.geometry
import os
from glob import glob


poni = r'C:\Users\kenneth1a\Documents\beamlineData\a311189_pdf_2/Si_0_15tilt.poni'
polarisation = 0.99
stdevs = 3
datadir = r'C:\Users\kenneth1a\Documents\beamlineData\a311189_pdf_2\01_RbCuF3_Q0p5\pdf\xrd/'
scale = 10**9

os.chdir(datadir)
if not os.path.exists(f'{datadir}/xye'):
    os.makedirs(f'{datadir}/xye/')

poni = pyFAI.load(poni)
geometry = pyFAI.geometry.Geometry(detector='pilatus2mcdte')
geometry.load(poni)


array2th = geometry.twoThetaArray()*180/np.pi
polarray  = geometry.polarization(factor = polarisation)

bins = 800
binsize = np.max(array2th)/bins
binarray = (array2th/binsize).astype(np.uint16)
cbfs = glob(f'*.cbf')
dataarray = fabio.open(cbfs[0]).data
binarray = np.where(dataarray < 0, -1, binarray)

mask = fabio.open(r'C:\Users\kenneth1a\Documents\beamlineData\Feb2023_gainMap/minMask.edf').data


for file in cbfs:
    print(file)
    dataarray = fabio.open(file).data
    xyefile = file.replace('.cbf','.xye')
    outfile = f'xye/{xyefile}'
    header = fabio.open(file).header["_array_data.header_contents"].split('\r\n# ')
    monitorCounts = int([line for line in header if 'Flux' in line][0].replace('Flux',''))
    
    array = np.empty(shape = dataarray.shape)
    for n in range(bins):
        array = np.where(binarray == n, dataarray, np.nan)/polarray
        stdev = np.nanstd(array)
        median = np.nanmedian(array)
        mask = np.where(array > median + stdevs*stdev, 1, mask)
    
    #fig,ax = plt.subplots(1,2)
    #ax[0].imshow(mask)
    #ax[1].imshow(dataarray,vmax = 2000)
    #plt.show()
    
    normArray = (dataarray/monitorCounts) * scale
    poni.integrate1d(normArray,mask = mask, filename = outfile, polarization_factor = polarisation,
                     unit = '2th_deg', correctSolidAngle = True, method = 'bbox', npt = 5000, 
                     error_model = 'poisson', safe = False)
    