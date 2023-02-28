#mask generator which goes through all directories in a given directory (unless 'average' is in the path)
import pyFAI
import fabio
import numpy as np
from glob import glob
import os


poni = r'C:\Users\kenneth1a\Documents\beamlineData\a311189_pdf_2/Si_0_15tilt.poni'
direc = r'C:\Users\kenneth1a\Documents\beamlineData\a311189_pdf_2/'
mask = r'C:\Users\kenneth1a\Documents\beamlineData\20220902_benchmark/pyfai-mask.edf'
gainFile = r'C:\Users\kenneth1a\Documents\maskGenerator/calculatedGainMap_48p6keV_kpm_filtered.edf'

os.chdir(direc)
mask = fabio.open(mask).data
poni = pyFAI.load(poni)



stdevs = 3
scale = 1e9
doMonitor = True
doGain = True
if doGain:
    gainArray = fabio.open(gainFile).data

for root, dirs, files in os.walk(direc):
    if 'average' in root:
        continue
    os.chdir(root)
    cbfs = glob('*.cbf')
    cbfs.sort()
    if len(cbfs) == 0:
        continue
    print(root)
    if doMonitor:
        monitorFile = glob('*.dat')[0]
        monitor = np.loadtxt(monitorFile,usecols = 2, skiprows = 1)
    subdir = f'xye_{stdevs}stdev/'
    subdirGain = f'xye_{stdevs}stdev_gainCorr/'
    if not os.path.exists(f'{root}/{subdir}/'):
        os.makedirs(f'{root}/{subdir}')
    if not os.path.exists(f'{root}/average/xye/'):
        os.makedirs(f'{root}/average/xye/')
    if not os.path.exists(f'{root}/{subdir}/average/'):
        os.makedirs(f'{root}/{subdir}/average/')
    if not os.path.exists(f'{root}/{subdirGain}/'):
        os.makedirs(f'{root}/{subdirGain}')
    if not os.path.exists(f'{root}/{subdirGain}/average/'):
        os.makedirs(f'{root}/{subdirGain}/average/')
    i1 = fabio.open(cbfs[0]).data
    dataset = np.empty(shape = (*i1.shape,len(cbfs)))
    for c,file in enumerate(cbfs):

        array = fabio.open(file).data
        if doMonitor:
            array = (array/monitor[c])*scale
        dataset[:,:,c] = array

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
    im.save(f'{root}/average/average.cbf')
    outfile = f'{root}/average/xye/average.xye'
    outfile_2d = outfile.replace('.xye','_pyfai.edf')
    mask_av = np.where(avim < 0, 1, 0)
    poni.integrate1d(data = avim, filename = outfile,mask = mask_av,polarization_factor = 0.99,unit = '2th_deg',
                correctSolidAngle = True, method = 'bbox',npt = 5000, error_model = 'poisson', safe = False)
    poni.integrate2d(data = avim, filename = outfile_2d,mask = mask_av,polarization_factor = 0.99,unit = '2th_deg',
                correctSolidAngle = True, method = 'bbox',npt_rad = 5000,npt_azim = 360, error_model = 'poisson', safe = False)
    if doGain:
        avimGain = avim/gainArray
        avimGain = np.where(gainArray <0, -1, avimGain)
        imGain = fabio.cbfimage.CbfImage(avimGain)
        imGain.save(f'{root}/average/average_gainCorrected.cbf')
    median = np.median(dataset, axis = 2)
    stdev = np.std(dataset,axis = 2)
    av1d = np.average(av1d,axis=1)
    eav = np.average(eav,axis = 1)
    np.savetxt(f'{subdir}average/{xyefile}',np.array([x,av1d,eav]).transpose())