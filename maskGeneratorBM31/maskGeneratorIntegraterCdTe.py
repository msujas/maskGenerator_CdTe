import pyFAI
import fabio
import numpy as np
from glob import glob
import os, re
import math

try:
    from . import  makeDataSet, makeMasks,  integrateAverage, integrateIndividual
except ImportError:
    from integrationFunctions import  makeDataSet, makeMasks,  integrateAverage, integrateIndividual

direc = r'Z:\visitor\a311217\bm31\20240129\pylatus\capillaries\LaB6_0-5_0-01\xrd/' # Directory of xrd files
dest = direc.replace(r'X:\users\a311207\20231204',r'Z:\visitor\a311207\bm31\20231204\pylatus')
poni  = r'Z:\visitor\a311217\bm31\20240129\pylatus/Si090_15tilt_MD.poni' # Poni file
mask  = r'Z:\visitor\a311217\bm31\20240129\pylatus/pdfMask.edf' # Mask file
gainFile = r'Z:\visitor\a311217\bm31\20240129\pylatus\gainmap\gainMap_thr30keV_filtered_kpm_2024-02-01.edf'
split = None

def run(direc,dest,poni,mask,gainFile, split = None, outdirav = 'average'):
    os.chdir(direc)
    if not os.path.exists(dest):
        os.makedirs(dest)

    mask = fabio.open(mask).data
    poni = pyFAI.load(poni)
    badFramesLog = f'{dest}/badFrames.txt'
    if os.path.isfile(badFramesLog):
        os.remove(badFramesLog)

    files = glob('*.cbf')
    if not files:
        return
    files.sort()
    if not split:
        splitVal = len(files)
    else:
        splitVal = split
    doMonitor = True
    scale = 10**9
    for i in range(math.ceil(len(files)/splitVal)):
        splitfiles = files[i*splitVal:(i+1)*splitVal]
        dataset, usedFiles = makeDataSet(splitfiles,badFramesLog,scale,doMonitor)
        print('\ncreating masks\n')
        nstdevs = 3
        maskdct = makeMasks(dataset, usedFiles, baseMask = mask, nstdevs = nstdevs)

        subdir = f'xye_{nstdevs}stdevs/'

        if not os.path.exists(f'{dest}/average/xye/'):
            os.makedirs(f'{dest}/average/xye/')

        print('\nmaking and integrating average images\n')
        integrateAverage(dataset, usedFiles, dest, poni, gainFile, maskdct, outdir=outdirav)
        print('\nintegrating individual images\n')
        integrateIndividual(dataset,usedFiles, dest, subdir, poni, maskdct, gainFile)
    
if __name__ == '__main__':
    run(direc=direc,dest=dest,poni=poni,mask=mask,gainFile=gainFile, split=split)
