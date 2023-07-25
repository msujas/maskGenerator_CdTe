import pyFAI
import fabio
import numpy as np
import matplotlib.pyplot as plt
from glob import glob
import os, re
from integrationFunctions import clearPyFAI_header, gainCorrection, bubbleHeader


direc = r'X:\staff\july2023\gunther\SiO2_4\xrd/' # Directory of xrd files
os.chdir(direc)

dest = direc.replace(r'X:\staff\july2023',r'C:\Users\kenneth1a\Documents\beamlineData\July2023')
#dest = direc.replace(r'X:\users\a311205',r'Z:\visitor\a311205\bm31\20230703\pylatus' )
if not os.path.exists(dest):
    os.makedirs(dest)

mask  = r'Z:\bm31\inhouse\july2023/pdf_baseMask_tilt.edf' # Mask file
mask = fabio.open(mask).data
poni  = r'Z:\bm31\inhouse\july2023/Si_15tilt_0p25579A.poni' # Poni file
poni = pyFAI.load(poni)
wavelength = poni.wavelength*10**10
gainFile = r'C:\Users\kenneth1a\Documents\beamlineData\July2023\gainmap\calculatedGainMap_48p6keV_filtered_kpm_2023-07-21.edf'

gainArray = fabio.open(gainFile).data


files = glob('*.cbf')
i1 = fabio.open(files[0]).data
dataset = np.empty(shape = (*i1.shape,len(files)))

doMonitor = True

#monitorfile = glob('*.dat')[0]
#monitor = np.loadtxt(monitorfile,usecols = 2, skiprows = 1)

scale = 10**9
for c,file in enumerate(files):
    
    array = fabio.open(file).data
    if doMonitor:
        fileheader = fabio.open(file).header["_array_data.header_contents"].split('\r\n#')
        monitor = int([item for item in fileheader if 'Flux' in item][0].replace('Flux',''))
        array = (array/monitor)*scale

    else:
        array = array*1000 #multiply by 1000 as 10^6 is common monitor value

    dataset[:,:,c] = array
    
average = np.average(dataset,axis=2)
median = np.median(dataset,axis=2)
stdev = np.std(dataset,axis = 2)

maskdct = {}
nstdevs = 3
for c in range(len(files)):
    print(files[c])
    array = dataset[:,:,c]
    vmax = np.percentile(np.where(np.isnan(array),0,array),99.5)
    maskdct[c] = np.where(array > median+nstdevs*stdev,1,mask)

subdir = f'xye_{nstdevs}stdevs/'
subdirGain = f'xye_{nstdevs}stdevs_gainCorrected/'
clearHeader = True

basefilename = os.path.basename(files[-1])
shortbasename = re.sub('_[0-9][0-9][0-9][0-9]p','',basefilename).replace('.cbf','')

if not os.path.exists(f'{dest}/average/xye/'):
    os.makedirs(f'{dest}/average/xye/')

dataset2 = np.empty(shape = dataset.shape)
for n in maskdct:
    dataset2[:,:,n] = np.where(maskdct[n] == 0, dataset[:,:,n], np.nan)
avim = np.nanmean(dataset2, axis = 2)
avim = np.where(np.isnan(avim), -2, avim)
im = fabio.cbfimage.CbfImage(avim)
im.save(f'{dest}/average/{shortbasename}_average.cbf')

avimGain = gainCorrection(avim,gainArray)
imGain = fabio.cbfimage.CbfImage(avimGain)
imGain.save(f'{dest}/average/{shortbasename}_average_gainCorrected.cbf')



outfile = f'{dest}/average/xye/{shortbasename}_average.xye'
outfile_2d = outfile.replace('.xye','_pyfai.edf')
mask_av = np.where(avim < 0, 1, 0)
mask_avGain = np.where(avimGain < 0, 1, 0)
poni.integrate1d(data = avim, filename = outfile,mask = mask_av,polarization_factor = 0.99,unit = '2th_deg',
                correctSolidAngle = True, method = 'bbox',npt = 5000, error_model = 'poisson', safe = False)
result = poni.integrate2d(data = avim, filename = outfile_2d,mask = mask_av,polarization_factor = 0.99,unit = '2th_deg',
                correctSolidAngle = True, method = 'bbox',npt_rad = 5000,npt_azim = 360, error_model = 'poisson', safe = False)
bubbleHeader(outfile_2d,*result[:3])

outfileGC = f'{dest}/average/xye/{shortbasename}_average_gainCorrected.xye'
poni.integrate1d(data = avimGain, filename = outfileGC,mask =mask_avGain,polarization_factor = 0.99,unit = '2th_deg',
                correctSolidAngle = True, method = 'bbox',npt = 5000, error_model = 'poisson', safe = False)
outfileGC_2d = outfileGC.replace('.xye','_pyfai.edf')
result = poni.integrate2d(data = avimGain, filename = outfileGC_2d,mask = mask_avGain,polarization_factor = 0.99,unit = '2th_deg',
                correctSolidAngle = True, method = 'bbox',npt_rad = 5000,npt_azim = 360, error_model = 'poisson', safe = False)
clearPyFAI_header(outfile)
clearPyFAI_header(outfileGC)
bubbleHeader(outfileGC_2d,*result[:3])


if not os.path.exists(f'{dest}/{subdir}/average/'):
    os.makedirs(f'{dest}/{subdir}/average/')
if not os.path.exists(f'{dest}/{subdirGain}/average/'):
    os.makedirs(f'{dest}/{subdirGain}/average/')  
for c,file in enumerate(files):
    print(file)
    xyefile = file.replace('.cbf','.xye')
    outputfile = f'{dest}/{subdir}/{xyefile}'
    x,y,e = poni.integrate1d(data = dataset[:,:,c], filename = outputfile,mask = maskdct[c],polarization_factor = 0.99,unit = '2th_deg',
                    correctSolidAngle = True, method = 'bbox',npt = 5000, error_model = 'poisson', safe = False)

    outputfileGC = f'{dest}/{subdirGain}/{xyefile}'
    arrayGC = dataset[:,:,c]/gainArray
    arrayGC = np.where(gainArray < 0, -1, arrayGC)
    xg,yg,eg = poni.integrate1d(data = arrayGC, filename = outputfileGC,mask = maskdct[c],polarization_factor = 0.99,unit = '2th_deg',
                    correctSolidAngle = True, method = 'bbox',npt = 5000, error_model = 'poisson', safe = False)
       
    if clearHeader:
        clearPyFAI_header(outputfile)
        clearPyFAI_header(outputfileGC)


    if c == 0:
        av1d = np.empty(shape = (len(y),len(files)))
        av1dg = np.empty(shape = (len(yg),len(files)))
        eav = np.empty(shape = (len(y),len(files)))
        eavg = np.empty(shape = (len(yg),len(files)))
    av1d[:,c] = y
    av1dg[:,c] = yg
    eav[:,c] = e
    eavg[:,c] = eg
av1d = np.average(av1d,axis=1)
av1dg = np.average(av1dg,axis=1)
eav = np.average(eav,axis = 1)
eavg = np.average(eavg,axis = 1)
np.savetxt(f'{dest}/{subdir}/average/{xyefile}',np.array([x,av1d,eav]).transpose())
np.savetxt(f'{dest}/{subdirGain}/average/{xyefile}',np.array([xg,av1dg,eavg]).transpose())
qav = av1d*4*np.pi*np.sin(x*np.pi/(180*2))/wavelength
qavg = av1dg*4*np.pi*np.sin(xg*np.pi/(180*2))/wavelength
np.savetxt(f'{dest}/{subdir}/average/qav.xy',np.array([x,qav]).transpose())
np.savetxt(f'{dest}/{subdirGain}/average/qav.xy',np.array([xg,qavg]).transpose())
