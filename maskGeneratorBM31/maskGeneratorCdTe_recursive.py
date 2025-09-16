#mask generator which goes through all directories in a given directory (unless 'average' is in the path)
import pyFAI
import fabio
import numpy as np
from glob import glob
import os
import math
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

def runSplit(root, stdevs, cbfs, poni, mask, gainFile, outfolder, splitval, badFramesLog, outdir, stop = False):
    subdir = f'xye_{stdevs}stdev/'
    if stop:
        return
    for i in range(math.ceil(len(cbfs)/splitval)):
        cbfstemp = cbfs[i*splitval:(i+1)*splitval]
        dataset, usedFiles = makeDataSet(cbfstemp,  badFramesLog, scale = scale, doMonitor = doMonitor)

        if np.any(dataset) == False:
            print(f'no monitor for {root}')
            continue
        print('\nmaking masks')
        masks = makeMasks(dataset = dataset,files =  usedFiles, baseMask = mask, nstdevs = 3)
        
        for n in range(5):
            try:
                print('\nmaking and integrating average image')
                integrateAverage(dataset, files = usedFiles, dest = outfolder, poni=poni, gainFile= gainFile, maskdct= masks, outdir = outdir)
                break
            except OSError as e:
                if n == 4:
                    raise OSError(e)
                time.sleep(1)
        print('\nintegrating individual images')
        integrateIndividual(dataset,files = usedFiles, dest = outfolder, subdir = subdir, avdir = avdir,  poni = poni, maskdct= masks, 
                            gainFile=gainFile)

def rundirSetup(root, basedir,dest,   fileList = None, split = None):
    if not fileList:
        fileList = []
    #if avdir in root or 'badFrames' in root or folderPattern not in root:
    #    return
    os.chdir(root)
    cbfs = glob('*.cbf')
    
    #if len(cbfs) == 0:
    #    return
    cbfs.sort()
    runDirec = False

    for cbf in cbfs:
        fullcbf = f'{root}/{cbf}'
        if not fullcbf in fileList:
            fileList.append(fullcbf)
            runDirec = True

    #if not runDirec:
    #    return
    
    #subdir = f'xye_{stdevs}stdev/'
    outfolder = root.replace(basedir,dest)
    badFramesLog = f'{outfolder}/badFrames.txt'
    if os.path.isfile(badFramesLog):
        os.remove(badFramesLog)
    splitval = split
    if not split:
        splitval = len(cbfs)
    return cbfs, outfolder, splitval, badFramesLog,  fileList, runDirec

def rundir(root, basedir,dest,poni,mask,gainFile, runningFull:bool,  fileList = None, split = None, 
           outdir = 'average', stop = False):

    cbfs, outfolder, splitval, badFramesLog, fileList, runDirec = rundirSetup(root=root, basedir=basedir,dest=dest, 
                                                                             fileList = fileList, split = split)
    
    if runDirec:
        runSplit(root, 3, cbfs, poni, mask, gainFile, outfolder, splitval, badFramesLog, outdir=outdir, stop = stop)
        runningFull = True
    return fileList, runningFull

def run(direc,dest,poni,mask,gainFile, folderPattern = '', fileList = None, split = None, outdir = 'average'):
    #if fileList == None:
    #    fileList = []
    #os.chdir(direc)
    mask = fabio.open(mask).data
    if isinstance(poni,str):
        poni = pyFAI.load(poni)
    
    runningFull = False
    for root, _dirs, _files in os.walk(direc):
        if outdir in root or 'badFrames' in root or folderPattern not in root:
            continue
        fileList, runningFull = rundir(root=root, direc=direc,dest=dest,poni=poni, mask=mask,gainFile=gainFile, 
                                       runningFull=runningFull, fileList=fileList, split=split, outdir=outdir)
    return fileList, runningFull
        
def runLoop(direc,dest,poni,mask,gainFile,folderPattern):
    fileList = []
    while True:
        fileList, runningFull = run(direc=direc,dest=dest,poni=poni,mask=mask, gainFile=gainFile, folderPattern=folderPattern, 
                                    fileList=fileList)
        time.sleep(1)
        if runningFull:
            print('looking for new files')


if __name__ == '__main__':
    run(direc=direc,dest=dest,poni=poni,mask=mask,gainFile=gainFile, folderPattern=folderPattern)    
     
