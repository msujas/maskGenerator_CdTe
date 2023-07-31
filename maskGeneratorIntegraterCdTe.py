import pyFAI
import fabio
import numpy as np
from glob import glob
import os, re
from integrationFunctions import  makeDataSet, makeMasks,  integrateAverage, integrateIndividual


direc = r'X:\staff\july2023\Julian\S22\xrd/' # Directory of xrd files
os.chdir(direc)

dest = direc.replace(r'X:\staff\july2023',r'C:\Users\kenneth1a\Documents\beamlineData\July2023')

if not os.path.exists(dest):
    os.makedirs(dest)

mask  = r'Z:\bm31\inhouse\july2023/pdf_baseMask_tilt2.edf' # Mask file
mask = fabio.open(mask).data
poni  = r'Z:\bm31\inhouse\july2023/Si_15tilt_0p25579A.poni' # Poni file
poni = pyFAI.load(poni)
wavelength = poni.wavelength*10**10
gainFile = r'C:\Users\kenneth1a\Documents\beamlineData\July2023\gainmap\calculatedGainMap_48p6keV_filtered_kpm_2023-07-21.edf'

gainArray = fabio.open(gainFile).data

badFramesLog = f'{dest}/badFrames.txt'
if os.path.isfile(badFramesLog):
    os.remove(badFramesLog)

files = glob('*.cbf')

doMonitor = True

#monitorfile = glob('*.dat')[0]
#monitor = np.loadtxt(monitorfile,usecols = 2, skiprows = 1)

scale = 10**9
dataset, usedFiles = makeDataSet(files,badFramesLog,scale,doMonitor)

    
average = np.average(dataset,axis=2)
median = np.median(dataset,axis=2)
stdev = np.std(dataset,axis = 2)
print('\ncreating masks\n')
nstdevs = 3
maskdct = makeMasks(dataset, usedFiles, baseMask = mask, nstdevs = nstdevs)

subdir = f'xye_{nstdevs}stdevs/'
subdirGain = f'xye_{nstdevs}stdevs_gainCorrected/'


if not os.path.exists(f'{dest}/average/xye/'):
    os.makedirs(f'{dest}/average/xye/')

print('\nmaking and integrating average images\n')
integrateAverage(dataset, usedFiles, poni, gainArray, maskdct, dest)
print('\nintegrating individual images\n')
integrateIndividual(dataset,usedFiles, dest, subdir, poni, maskdct, gainArray)

