#mask generator which goes through all directories in a given directory (unless 'average' is in the path)
import pyFAI
import fabio
import numpy as np
from glob import glob
import os
import math
from pathlib import Path
import time
from maskGeneratorBM31 import  Integrator

class RecursiveIntegrator():
    def __init__(self, direc, stdevs, poni, mask, gainfile,   saveindividual=True, splitval=0, polF = 0.85, avdir = 'average',
                 folderpattern=''):
        self.direc = direc
        self.stdevs = stdevs
        self.poni=poni
        self.mask=mask
        self.gainfile=gainfile
        self.splitval=splitval
        self.polF=polF
        self.avdir = avdir
        self.filelist= []
        self.saveindividual=saveindividual
        self.folderpattern=folderpattern
        self.running = True
        

    def runSplit(self, root:str, scale=10**9, doMonitor=True):
        subdir = f'xye_{self.stdevs}stdev/'
        if not self.running:
            return
        if not self.splitval:
            splitval = len(self.cbfs)
        else:
            splitval = self.splitval
        for i in range(math.ceil(len(self.cbfs)/splitval)):
            cbfstemp = self.cbfs[i*splitval:(i+1)*splitval]
            self.integrator.files = cbfstemp
            self.integrator.makeDataSet(scale = scale, doMonitor = doMonitor)

            if np.any(self.integrator.dataset) == False:
                print(f'no monitor for {root}')
                continue
            print('\nmaking masks')
            self.integrator.makeMasks()
            
            for n in range(5):
                try:
                    print('\nmaking and integrating average image')
                    self.integrator.integrateAverage()
                    break
                except OSError as e: #if file not on network yet
                    if n == 4:
                        raise OSError(e)
                    time.sleep(1)
            if self.saveindividual:
                print('\nintegrating individual images')
                self.integrator.integrateIndividual(subdir = subdir)

    def rundirSetup(self,root):
        self.cbfs = glob(f'{root}/*.cbf')
        self.cbfs.sort()
        runDirec = False

        for cbf in self.cbfs:
            if not cbf in self.filelist:
                self.filelist.append(cbf)
                runDirec = True

        return runDirec

    def rundir(self,root):
        runDirec = self.rundirSetup(root=root)
        if runDirec:
            self.runSplit(root)
            self.runningFull = True

    def run(self):
        self.runningFull = False
        for root, _dirs, _files in os.walk(self.direc):
            self.integrator = Integrator(root, self.mask, self.poni, self.gainfile, self.polF, avdir=self.avdir)
            if self.avdir in root or 'badFrames' in root or self.folderpattern not in root:
                continue
            self.rundir(root)

            
    def runLoop(self):
        while self.running:
            self.run()
            time.sleep(1)
            if self.runningFull:
                print('looking for new files')
