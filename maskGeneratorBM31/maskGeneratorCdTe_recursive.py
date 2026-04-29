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

class RecursiveIntegrator():
    def __init__(self, root, stdevs, poni, mask, gainfile, dest, splitval=0):
        self.root = root
        self.stdevs = stdevs
        self.poni=poni
        self.mask=mask
        self.gainfile=gainfile
        self.outfolder=dest
        self.splitval=splitval
        self.cbfs = glob(f'{root}/*.cbf')
        
        if not self.splitval:
            self.splitval = len(self.cbfs)
        self.filelist= []


        

    def runSplit(self, stdevs,  poni, mask, gainFile, outfolder, badFramesLog, outdir, stop = False, individual = True,
                polF = 0.85):
        subdir = f'xye_{stdevs}stdev/'
        if stop:
            return
        for i in range(math.ceil(len(self.cbfs)/self.splitval)):
            cbfstemp = self.cbfs[i*self.splitval:(i+1)*self.splitval]
            dataset, usedFiles = makeDataSet(cbfstemp,  badFramesLog, scale = scale, doMonitor = doMonitor)

            if np.any(dataset) == False:
                print(f'no monitor for {self.root}')
                continue
            print('\nmaking masks')
            masks = makeMasks(dataset = dataset,files =  usedFiles, baseMask = mask, nstdevs = 3)
            
            for n in range(5):
                try:
                    print('\nmaking and integrating average image')
                    integrateAverage(dataset, files = usedFiles, dest = outfolder, poni=poni, gainFile= gainFile, maskdct= masks, outdir = outdir, polF=polF)
                    break
                except OSError as e:
                    if n == 4:
                        raise OSError(e)
                    time.sleep(1)
            if individual:
                print('\nintegrating individual images')
                integrateIndividual(dataset,files = usedFiles, dest = outfolder, subdir = subdir, avdir = avdir,  poni = poni, maskdct= masks, 
                                gainFile=gainFile, polF=polF)

    def rundirSetup(self,root):

        #if avdir in root or 'badFrames' in root or folderPattern not in root:
        #    return
        cbfs = glob(f'{root}/*.cbf')
        
        #if len(cbfs) == 0:
        #    return
        cbfs.sort()
        runDirec = False

        for cbf in cbfs:
            if not cbf in self.filelist:
                self.filelist.append(cbf)
                runDirec = True

        #if not runDirec:
        #    return
        
        #subdir = f'xye_{stdevs}stdev/'

        return cbfs, runDirec

    def rundir(self,root, basedir,dest,poni,mask,gainFile, runningFull:bool,  split = None, 
            outdir = 'average', stop = False, individual = True, polF = 0.99):

        cbfs,   runDirec = self.rundirSetup(root=root, basedir=basedir,dest=dest, 
                                            split = split)
        
        if runDirec:
            self.runSplit(root, 3, cbfs, poni, mask, gainFile,  outdir=outdir, stop = stop, individual=individual,polF=polF)
            runningFull = True
        return runningFull

    def run(self,direc,dest,poni,mask,gainFile, folderPattern = '', fileList = None, split = None, outdir = 'average'):
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
            fileList, runningFull = self.rundir(root=root, direc=direc,dest=dest,poni=poni, mask=mask,gainFile=gainFile, 
                                        runningFull=runningFull, fileList=fileList, split=split, outdir=outdir)
        return fileList, runningFull
            
    def runLoop(self,direc,dest,poni,mask,gainFile,folderPattern):
        while True:
            fileList, runningFull = self.run(direc=direc,dest=dest,poni=poni,mask=mask, gainFile=gainFile, folderPattern=folderPattern)
            time.sleep(1)
            if runningFull:
                print('looking for new files')


if __name__ == '__main__':
    run(direc=direc,dest=dest,poni=poni,mask=mask,gainFile=gainFile, folderPattern=folderPattern)    
     
