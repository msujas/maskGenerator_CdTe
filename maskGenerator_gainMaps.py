import pyFAI
import fabio
import numpy as np
from glob import glob
import os, re
import maskGeneratorBM31 as ifuns


direc =  r'C:\Users\kenneth1a\Documents\beamlineData\March2023_gainMap\GainMeasure\pos7/xrd/' # Directory of xrd files
mask  = fr'C:\Users\kenneth1a\Documents\beamlineData\Feb2023_gainMap\minMask.edf' # Mask file
gainFile = r'C:\Users\kenneth1a\Documents\beamlineData\Feb2023_gainMap\glass2\minMask/calculatedGainMap_48p6keV.edf'
ponidir = r'C:\Users\kenneth1a\Documents\beamlineData\March2023_gainMap/'
ponis = glob(f'{ponidir}/*.poni')

pos = [x for x in re.split('/|\\\\',direc) if 'pos' in x][0]
poni = [poni for poni in ponis if pos in poni][0]
def run(direc,mask,gainFile,poni,avdir = 'average', dest = None,doMonitor = True):
    '''
    direc - directory containing cbf files
    mask - base mask for file for integrateion
    gainFile - file with gain values of pixels
    avdir - destination directory for averaged images and integrations of them
    '''
    
    print(poni)
    os.chdir(direc)
    if dest == None:
        dest = direc
    if not os.path.exists(dest):
        os.makedirs(dest)


    mask = fabio.open(mask).data

    poni = pyFAI.load(poni)

    files = glob('*.cbf')
    files.sort()
    i1 = fabio.open(files[0]).data
    dataset = np.empty(shape = (*i1.shape,len(files)))

    scale = 10**9
    badFramesLog = f'{dest}/badFrames.log'
    dataset,usedfiles = ifuns.makeDataSet(files,badFramesLog, scale, doMonitor=doMonitor)

    nstdevs = 3
    maskdct = ifuns.makeMasks(dataset, usedfiles,mask,nstdevs)

    subdir = f'xye_{nstdevs}stdevs/'

    ifuns.integrateAverage(dataset, usedfiles, dest, poni, gainFile, maskdct, unit = '2th_deg', npt = 5000, nptA = 360, polF = 0.99)
    ifuns.integrateIndividual(dataset,usedfiles, dest, subdir, poni, maskdct, gainFile, avdir = 'average', unit = '2th_deg', npt = 5000, polF = 0.99)

if __name__ == '__main__':
    run(direc,mask,gainFile, poni)