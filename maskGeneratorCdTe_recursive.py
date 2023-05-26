#mask generator which goes through all directories in a given directory (unless 'average' is in the path)
import pyFAI
import fabio
import numpy as np
from glob import glob
import os, re
import matplotlib.pyplot as plt
from integrationFunctions import clearPyFAI_header, gainCorrection, bubbleHeader


poni = r'C:\Users\kenneth1a\Documents\beamlineData\May2023\Si00_tilt/pos05.poni'
direc = r'C:\Users\kenneth1a\Documents\beamlineData\May2023\measurements'
mask = r'C:\Users\kenneth1a\Documents\beamlineData\May2023\Si00_tilt\pos05_mask.edf'
gainFile = r'C:\Users\kenneth1a\Documents\beamlineData\May2023\compare/afterRestart_solidAngle.edf'
#gainFile = r'C:\Users\kenneth1a\Documents\beamlineData\Feb2023_gainMap/calculatedGainMap_48p6keV_2023-02-21_2.edf'
#gainFile = r'C:\Users\kenneth1a\Documents\beamlineData\Feb2023_gainMap/calculatedGainMap_48p6keV_filtered_kpm_2023-04-13.edf'
os.chdir(direc)
mask = fabio.open(mask).data
poni = pyFAI.load(poni)


avdir = 'average'
stdevs = 3
scale = 1e9
doMonitor = True
doGain = True
if doGain:
    gainArray = fabio.open(gainFile).data

for root, dirs, files in os.walk(direc):
    if avdir in root:
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
    subdirGain = f'xye_{stdevs}stdev_gainCorr/'
    if not os.path.exists(f'{root}/{subdir}/'):
        os.makedirs(f'{root}/{subdir}')
    if not os.path.exists(f'{root}/{avdir}/xye/'):
        os.makedirs(f'{root}/{avdir}/xye/')
    if not os.path.exists(f'{root}/{subdir}/{avdir}/'):
        os.makedirs(f'{root}/{subdir}/{avdir}/')
    if not os.path.exists(f'{root}/{subdirGain}/'):
        os.makedirs(f'{root}/{subdirGain}')
    if not os.path.exists(f'{root}/{subdirGain}/{avdir}/'):
        os.makedirs(f'{root}/{subdirGain}/{avdir}/')
    i1 = fabio.open(cbfs[0]).data
    dataset = np.empty(shape = (*i1.shape,len(cbfs)))
    for c,file in enumerate(cbfs):

        array = fabio.open(file).data
        if doMonitor:
            try:
                #monitor = monitorList[c] #if monitor in separate file instead of header
                fileheader = fabio.open(file).header["_array_data.header_contents"].split('\r\n#')
                monitor = int([item for item in fileheader if 'Flux' in item][0].replace('Flux',''))
                
            except:
                continue
            array = (array/monitor)*scale  
        dataset[:,:,c] = array
    if np.any(dataset) == False:
        print(f'no monitor for {root}')
        continue
    median = np.median(dataset,axis = 2)
    stdev = np.std(dataset,axis = 2)

    masks = {}
    for c,file in enumerate(cbfs):
        array = dataset[:,:,c]
        masks[c] = np.where(array > median+stdevs*stdev,1,mask)

        xyefile = file.replace('.cbf','.xye')
        outputfile = f'{root}/{subdir}/{xyefile}'
        x,y,e = poni.integrate1d(data = dataset[:,:,c], filename = outputfile,mask = masks[c],polarization_factor = 0.99,unit = '2th_deg',
                    correctSolidAngle = True, method = 'bbox',npt = 5000, error_model = 'poisson',safe = False)
        print(outputfile)
        if doGain:
            arrayGC = dataset[:,:,c]/gainArray
            arrayGC = np.where(gainArray < 0, -1, arrayGC)
            outputfileGC = f'{root}/{subdirGain}/{xyefile}'
            xg,yg,eg = poni.integrate1d(data = arrayGC, filename = outputfileGC,mask = masks[c],polarization_factor = 0.99,unit = '2th_deg',
            correctSolidAngle = True, method = 'bbox',npt = 5000, error_model = 'poisson',safe = False)
        if c == 0:
            av1d = np.empty(shape = (len(y),len(files)))
            eav = np.empty(shape = (len(y),len(files)))
            av1d2 = np.empty(shape = (len(y),len(files)))
            eav2 = np.empty(shape = (len(y),len(files)))
        av1d[:,c] = y
        eav[:,c] = e
    dataset2 = np.empty(shape = dataset.shape)
    for n in masks:
        dataset2[:,:,n] = np.where(masks[n] == 0, dataset[:,:,n], np.nan)
    avim = np.nanmean(dataset2, axis = 2)
    avim = np.where(np.isnan(avim), -2, avim)
    im = fabio.cbfimage.CbfImage(avim)
    im.save(f'{root}/{avdir}/average.cbf')

    mask_av = np.where(avim < 0, 1, 0)
    basefilename = os.path.basename(cbfs[-1])
    shortbasename = re.sub('_[0-9][0-9][0-9][0-9]p','',basefilename).replace('.cbf','')
    outfile = f'{root}/{avdir}/xye/{shortbasename}_average.xye'
    outfile_2d = outfile.replace('.xye','_pyfai.edf')
    poni.integrate1d(data = avim, filename = outfile,mask = mask_av,polarization_factor = 0.99,unit = '2th_deg',
                correctSolidAngle = True, method = 'bbox',npt = 5000, error_model = 'poisson', safe = False)
    array2d, tth, eta, e = poni.integrate2d(data = avim, filename = outfile_2d,mask = mask_av,polarization_factor = 0.99,unit = '2th_deg',
                correctSolidAngle = True, method = 'bbox',npt_rad = 5000,npt_azim = 360, error_model = 'poisson', safe = False)
    clearPyFAI_header(outfile)
    bubbleHeader(outfile_2d, array2d, tth, eta)
    if doGain:
        avimGain = gainCorrection(avim,gainArray)
        imGain = fabio.cbfimage.CbfImage(avimGain)
        imGain.save(f'{root}/{avdir}/average_gainCorrected.cbf')
        outfileGC = f'{root}/{avdir}/xye/{shortbasename}_average_gainCorrected.xye'
        outfile_2dGC = outfileGC.replace('.xye','.edf')
        poni.integrate1d(data = avimGain, filename = outfileGC,mask = mask_av,polarization_factor = 0.99,unit = '2th_deg',
            correctSolidAngle = True, method = 'bbox',npt = 5000, error_model = 'poisson', safe = False)
        array2d, tth, eta, e = poni.integrate2d(data = avimGain, filename = outfile_2dGC,mask = mask_av,polarization_factor = 0.99,unit = '2th_deg',
                correctSolidAngle = True, method = 'bbox',npt_rad = 5000,npt_azim = 360, error_model = 'poisson', safe = False)
        bubbleHeader(outfile_2dGC, array2d, tth, eta)
        clearPyFAI_header(outfileGC)      
    median = np.median(dataset, axis = 2)
    stdev = np.std(dataset,axis = 2)
    av1d = np.average(av1d,axis=1)
    eav = np.average(eav,axis = 1)
    np.savetxt(f'{subdir}/{avdir}/{xyefile}',np.array([x,av1d,eav]).transpose())