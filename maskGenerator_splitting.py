#mask generator which splits large datasets up, creating multiple averages in addition to the individual patterns.
#useful for datasets with 100s or 1000s of images which would require too much storage to do the mask generation 
#for everything at once

import pyFAI
import fabio
import numpy as np
from glob import glob
import os
from cryio.cbfimage import CbfHeader
from maskGeneratorBM31 import Integrator


#direc = os.getcwd() # Current Directory
direc = r'X:\users\a311207\20231204\Basset_lfp2\xrd/' # Directory of xrd files
dest = direc.replace(r'X:\users\a311207\20231204',r'Z:\visitor\a311207\bm31\20231204\pylatus')

poni = r'Z:\visitor\a311207\bm31\20231204\pylatus/Si000_15tilt.poni'
mask = r'Z:\visitor\a311207\bm31\20231204\pylatus\pdf_baseMask_tilt.edf'
gainFile = r'Z:\visitor\a311207\bm31\20231204\pylatus/calculatedGainMap_48p6keV_filtered_kpm_2023-12-04.edf'
pfactor = 0.85

split = 10 #change this depending on how many you want to average

def run(direc, dest,poni,mask,gainFile,split = 20, pfactor = 0.85):
    i = Integrator(direc, mask, poni, gainFile, polfactor=pfactor, dest=dest)
    i.runsplit(split=split)
        
if __name__ == '__main__':
    run(direc,dest,poni,mask,gainFile,split)