#mask generator which goes through all directories in a given directory (unless 'average' is in the path)
import pyFAI
import fabio
import numpy as np
from glob import glob
import os
from pathlib import Path
import time
from maskGeneratorBM31 import  makeDataSet, integrateAverage, integrateIndividual, makeMasks


direc = r'X:\users\a311207\20231204'
homedir = str(Path.home())

dest = direc.replace(r'X:\users\a311207\20231204',fr'Z:\visitor\a311207\bm31\20231204\pylatus')
poni = r'X:\users\a311207\20231204/Si_0_15tilt.poni'
mask = r'X:\users\a311207\20231204\dtx0_dtr15_baseMask_lines.edf'
gainFile =  r'X:\users\a311207\20231204/calculatedGainMap_48p6keV_kpm_filtered.edf'
avdir = 'average'
stdevs = 3
scale = 1e9
doMonitor = True
folderPattern = 'pdf' #pattern to search for somewhere in the directory name

def run(direc,dest,poni,mask,gainFile, folderPattern = '', fileList = None):
    if fileList == None:
        fileList = []
    os.chdir(direc)
    mask = fabio.open(mask).data
    poni = pyFAI.load(poni)
    runningFull = False
    for root, dirs, files in os.walk(direc):
        if avdir in root or 'badFrames' in root or folderPattern not in root:
            continue
        os.chdir(root)
        cbfs = glob('*.cbf')
        cbfs.sort()
        if len(cbfs) == 0:
            continue
        
        runDirec = False
        for cbf in cbfs:
            fullcbf = f'{root}/{cbf}'
            if not fullcbf in fileList:
                fileList.append(fullcbf)
                runDirec = True
                runningFull = True
        if not runDirec:
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
        for n in range(5):
            try:
                integrateAverage(dataset, files = usedFiles, dest = outfolder, poni=poni, gainFile= gainFile, maskdct= masks)
                print('\nintegrating individual images')
            except OSError as e:
                if n == 4:
                    raise OSError(e)
                time.sleep(1)
        integrateIndividual(dataset,files = usedFiles, dest = outfolder, subdir = subdir, avdir = avdir,  poni = poni, maskdct= masks, 
                            gainFile=gainFile)
    return fileList, runningFull
        
def runLoop(direc,dest,poni,mask,gainFile,folderPattern):
    fileList = []
    while True:
        fileList, runningFull = run(direc,dest,poni,mask, gainFile, folderPattern, fileList)
        time.sleep(1)
        if runningFull:
            print('looking for new files')


if __name__ == '__main__':
    run(direc=direc,dest=dest,poni=poni,mask=mask,gainFile=gainFile, folderPattern=folderPattern)    
     
