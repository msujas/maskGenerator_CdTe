#mask generator which goes through all directories in a given directory (unless 'average' is in the path)
import pyFAI
import fabio
import numpy as np
from glob import glob
import os
from pathlib import Path
from maskGeneratorBM31 import  makeDataSet, integrateAverage, integrateIndividual, makeMasks

direc = r'X:\users\a311207\20231204\capillaries'
homedir = str(Path.home())

dest = direc.replace(r'X:\users\a311207\20231204',fr'Z:\visitor\a311207\bm31\20231204\pylatus')
poni = r'Z:\visitor\a311207\bm31\20231204\pylatus/Si000_15tilt.poni'
mask = r'Z:\visitor\a311207\bm31\20231204\pylatus\pdf_baseMask_tilt.edf'
gainFile =  fr'Z:\visitor\a311207\bm31\20231204\pylatus/calculatedGainMap_48p6keV_filtered_kpm_2023-12-04.edf'
avdir = 'average'
stdevs = 3
scale = 1e9
doMonitor = True



def run(direc,dest,poni,mask,gainFile):
    os.chdir(direc)
    mask = fabio.open(mask).data
    poni = pyFAI.load(poni)


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
        print('\nmaking masks')
        masks = makeMasks(dataset = dataset,files =  usedFiles, baseMask = mask, nstdevs = 3)
        print('\nmaking and integrating average image')
        integrateAverage(dataset, files = usedFiles, dest = outfolder, poni=poni, gainFile= gainFile, maskdct= masks)
        print('\nintegrating individual images')
        integrateIndividual(dataset,files = usedFiles, dest = outfolder, subdir = subdir, avdir = avdir,  poni = poni, maskdct= masks, 
                            gainFile=gainFile)
if __name__ == '__main__':
    run(direc=direc,dest=dest,poni=poni,mask=mask,gainFile=gainFile)    
     
