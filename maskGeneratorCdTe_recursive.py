#mask generator which goes through all directories in a given directory (unless 'average' is in the path)
import pyFAI
import fabio
import numpy as np
from glob import glob
import os
from integrationFunctions import  makeDataSet, integrateAverage, integrateIndividual, makeMasks

direc = r'X:\staff\july2023\Julian'
dest = r'C:\Users\kenneth1a\Documents\beamlineData\July2023\Julian'

poni = r'C:\Users\kenneth1a\Documents\beamlineData\July2023/Si_15tilt_0p25579A.poni'
mask = r'C:\Users\kenneth1a\Documents\beamlineData\July2023\pdf_baseMask_tilt2.edf'
gainFile = r'C:\Users\kenneth1a\Documents\beamlineData\July2023\gainmap/calculatedGainMap_48p6keV_filtered_kpm_2023-07-21.edf'


os.chdir(direc)
mask = fabio.open(mask).data
poni = pyFAI.load(poni)


avdir = 'average'
stdevs = 3
scale = 1e9
doMonitor = True
doGain = True
if doGain:
    gainArray = fabio.open(gainFile).data

for root, dirs, files in os.walk(direc):
    if avdir in root or 'badFrames' in root:
        continue
    os.chdir(root)
    cbfs = glob('*.cbf')
    cbfs.sort()
    if len(cbfs) == 0:
        continue
    print(root)
    #if doMonitor: #for using monitor log files instead of file header
        #monitorFile = glob('*.dat')[0]
        #monitorList = np.loadtxt(monitorFile,usecols = 2, skiprows = 1)

    subdir = f'xye_{stdevs}stdev/'
    outfolder = root.replace(direc,dest)
    badFramesLog = f'{outfolder}/badFrames.txt'
    if os.path.isfile(badFramesLog):
        os.remove(badFramesLog)


    dataset, usedFiles = makeDataSet(cbfs,  badFramesLog, scale = scale, doMonitor = doMonitor)

    if np.any(dataset) == False:
        print(f'no monitor for {root}')
        continue
    median = np.median(dataset,axis = 2)
    stdev = np.std(dataset,axis = 2)

    masks = makeMasks(dataset = dataset,files =  usedFiles, baseMask = mask, nstdevs = 3)
    integrateAverage(dataset, files = usedFiles, dest = outfolder, poni=poni, gainArray= gainArray, maskdct= masks)
    integrateIndividual(dataset,files = usedFiles, dest = outfolder, subdir = subdir, avdir = avdir,  poni = poni, maskdct= masks, 
                        gainArray=gainArray)
    
     
