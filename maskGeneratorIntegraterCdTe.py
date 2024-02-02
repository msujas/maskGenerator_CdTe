import pyFAI
import fabio
import numpy as np
from glob import glob
import os, re
from integrationFunctions import  makeDataSet, makeMasks,  integrateAverage, integrateIndividual


direc = r'X:\users\a311207\20231204\air2\xrd/' # Directory of xrd files
dest = direc.replace(r'X:\users\a311207\20231204',r'Z:\visitor\a311207\bm31\20231204\pylatus')
poni  = r'Z:\visitor\a311207\bm31\20231204\pylatus/Si000_15tilt.poni' # Poni file
mask  = r'Z:\visitor\a311207\bm31\20231204\pylatus/pdf_baseMask_tilt.edf' # Mask file
gainFile = r'Z:\visitor\a311207\bm31\20231204\pylatus\calculatedGainMap_48p6keV_filtered_kpm_2023-12-11.edf'

def run(direc,dest,poni,mask,gainFile):
    os.chdir(direc)
    if not os.path.exists(dest):
        os.makedirs(dest)

    mask = fabio.open(mask).data
    poni = pyFAI.load(poni)
    badFramesLog = f'{dest}/badFrames.txt'
    if os.path.isfile(badFramesLog):
        os.remove(badFramesLog)

    files = glob('*.cbf')
    files.sort()
    doMonitor = True
    scale = 10**9
    dataset, usedFiles = makeDataSet(files,badFramesLog,scale,doMonitor)
    print('\ncreating masks\n')
    nstdevs = 3
    maskdct = makeMasks(dataset, usedFiles, baseMask = mask, nstdevs = nstdevs)

    subdir = f'xye_{nstdevs}stdevs/'


    if not os.path.exists(f'{dest}/average/xye/'):
        os.makedirs(f'{dest}/average/xye/')

    print('\nmaking and integrating average images\n')
    integrateAverage(dataset, usedFiles, dest, poni, gainFile, maskdct)
    print('\nintegrating individual images\n')
    integrateIndividual(dataset,usedFiles, dest, subdir, poni, maskdct, gainFile)
    
if __name__ == '__main__':
    run(direc=direc,dest=dest,poni=poni,mask=mask,gainFile=gainFile)
