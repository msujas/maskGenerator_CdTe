#mask generator which splits large datasets up, creating multiple averages in addition to the individual patterns.
#useful for datasets with 100s or 1000s of images which would require too much storage to do the mask generation 
#for everything at once

import pyFAI
import fabio
import numpy as np
from glob import glob
import os
from integrationFunctions import makeDataSet, makeMasks, integrateAverage, integrateIndividual


#direc = os.getcwd() # Current Directory
direc = r'X:\users\a311207\20230412\Basset_lfp2\xrd/' # Directory of xrd files
dest = direc.replace(r'X:\users\a311207\20230412',r'Z:\visitor\a311207\bm31\20231204\pylatus')

poni = r'Z:\visitor\a311207\bm31\20231204\pylatus/Si000_15tilt.poni'
mask = r'Z:\visitor\a311207\bm31\20231204\pylatus\pdf_baseMask_tilt.edf'
gainFile = r'Z:\visitor\a311207\bm31\20231204\pylatus/calculatedGainMap_48p6keV_filtered_kpm_2023-12-04.edf'


averaging = 10 #change this depending on how many you want to average

def run(direc, dest,poni,mask,gainFile,averaging = 20,doMonitor = True):
    os.chdir(direc)
    
    if not os.path.exists(dest):
        os.makedirs(dest)



    poni = pyFAI.load(poni)
    wavelength = poni.wavelength*10**10
   
    nstdevs = 3 #change this to make pixel masking stricter or more lenient
    scale = 10**9 #scaling monitor normalised data 10^9 should be good but can be adjusted

    files = glob('*.cbf')
    files.sort()
    filesplit = []
    n = -1

    for c,f in enumerate(files):
        if c%averaging == 0:
            filesplit.append([])
            n += 1
        filesplit[n].append(f)
            
    gainArray = fabio.open(gainFile).data
    mask = fabio.open(mask).data
    subdir = f'xye_{nstdevs}stdevs_{averaging}/'

    badFramesLog = f'{dest}/badFrames.log'
    for i,files in enumerate(filesplit):
        
        dataset, usedFiles = makeDataSet(files, badFramesLog, scale = scale, doMonitor = True)
            
        maskdct = makeMasks(dataset, usedFiles, mask, nstdevs = 3, plot = False)
        
        integrateAverage(dataset, usedFiles, dest, poni, gainArray, maskdct, unit = '2th_deg', npt = 5000, nptA = 360, polF = 0.99, 
                         fileAppend = f'_{i:04d}')
        
        integrateIndividual(dataset,usedFiles, dest, subdir, poni, maskdct, gainArray, avdir = 'average', unit = '2th_deg', 
                            npt = 5000, polF = 0.99)
        
if __name__ == '__main__':
    run(direc,dest,poni,mask,gainFile,averaging)