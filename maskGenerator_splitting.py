#mask generator which splits large datasets up, creating multiple averages in addition to the individual patterns.
#useful for datasets with 100s or 1000s of images which would require too much storage to do the mask generation 
#for everything at once

import pyFAI
import fabio
import numpy as np
from glob import glob
import os, re
from integrationFunctions import clearPyFAI_header, gainCorrection, bubbleHeader


#direc = os.getcwd() # Current Directory
direc = r'C:\Users\kenneth1a\Documents\beamlineData\Feb2023_gainMap\glass2\pos1/' # Directory of xrd files
poni  = fr'C:\Users\kenneth1a\Documents\beamlineData\Feb2023_gainMap\Si\xrd/Si_700_dty_179.poni' # Poni file
mask  = r'C:\Users\kenneth1a\Documents\beamlineData\Feb2023_gainMap\minMask.edf' # Mask file
gainFile = r'C:\Users\kenneth1a\Documents\maskGenerator/calculatedGainMap_48p6keV_kpm_filtered.edf'
averaging = 20 #change this depending on how many you want to average
def run(direc,poni,mask,gainFile,averaging = 20,doMonitor = True):
    os.chdir(direc)


    dest = direc
    if not os.path.exists(dest):
        os.makedirs(dest)

    mask = fabio.open(mask).data

    poni = pyFAI.load(poni)
    wavelength = poni.wavelength*10**10
   
    nstdevs = 3 #change this to make pixel masking stricter or more lenient
    scale = 10**9 #scaling monitor normalised data 10^9 should be good but can be adjusted

    files = glob('*.cbf')
    files.sort()
    filesplit = []
    n = -1

    for c,f in enumerate(files):
        if c%averaging == 0:
            filesplit.append([])
            n += 1
        filesplit[n].append(f)
            
    i1 = fabio.open(files[0]).data




    gainArray = fabio.open(gainFile).data


    avdir = f'average{averaging}/'
    subdir = f'xye_{nstdevs}stdevs_{averaging}/'
    subdirGain = f'xye_{nstdevs}stdevs_gainCorr_{averaging}/'
    clearHeader = True

    if not os.path.exists(f'{dest}/{avdir}/xye'):
        os.makedirs(f'{dest}/{avdir}/xye')
    if not os.path.exists(f'{dest}/{subdir}/{avdir}/'):
        os.makedirs(f'{dest}/{subdir}/{avdir}/')
    if not os.path.exists(f'{dest}/{subdirGain}/{avdir}/'):
        os.makedirs(f'{dest}/{subdirGain}/{avdir}/') 

    for i,files in enumerate(filesplit):
        dataset = np.empty(shape = (*i1.shape,len(files)))
        for c,file in enumerate(files):
            array = fabio.open(file).data
            if doMonitor:
                fileheader = fabio.open(file).header["_array_data.header_contents"].split('\r\n#')
                fluxline = [line for line in fileheader if 'Flux' in line][0]
                monitor = int(fluxline.replace('Flux',''))
                array = (array/monitor)*scale
            dataset[:,:,c] = array
            

        average = np.average(dataset,axis=2)
        median = np.median(dataset,axis=2)

        stdev = np.std(dataset,axis = 2)


        maskdct = {}

        for c in range(len(files)):
            print(files[c])
            array = dataset[:,:,c]
            maskdct[c] = np.where(array > median+nstdevs*stdev,1,mask)
        

        
        dataset2 = np.empty(shape = dataset.shape)
        for n in maskdct:
            dataset2[:,:,n] = np.where(maskdct[n] == 0, dataset[:,:,n], np.nan)
        avim = np.nanmean(dataset2, axis = 2)
        avim = np.where(np.isnan(avim), -2, avim)
        im = fabio.cbfimage.CbfImage(avim)
        im.save(f'{dest}/{avdir}/average_{i}.cbf')
        

        avimGain = gainCorrection(avim,gainArray)
        imGain = fabio.cbfimage.CbfImage(avimGain)
        imGain.save(f'{dest}/{avdir}/average_gainCorrected_{i}.cbf')
        
        basefilename = os.path.basename(files[-1])
        shortbasename = re.sub('_[0-9][0-9][0-9][0-9]p','',basefilename).replace('.cbf','')
        
        outfile = f'{dest}/{avdir}/xye/{shortbasename}_average_{i}.xye'
        outfile2d = outfile.replace('.xye','_pyfai.edf')
        mask_av = np.where(avim <0, 1, 0)

        poni.integrate1d(data = avim, filename = outfile,mask = mask_av,polarization_factor = 0.99,unit = '2th_deg',
                        correctSolidAngle = True, method = 'bbox',npt = 5000, error_model = 'poisson', safe = False)
        pattern2d, tth, eta, e = poni.integrate2d(data = avim, filename = outfile2d,mask = mask_av,polarization_factor = 0.99,unit = '2th_deg',
                        correctSolidAngle = True, method = 'bbox',npt_rad = 5000, npt_azim = 360, error_model = 'poisson', safe = False)
        bubbleHeader(outfile2d, pattern2d, tth, eta)
        outfileGain = f'{dest}/{avdir}/xye/average_gainCorrected_{i}.xye'
        outfileGain2d = outfileGain.replace('.xye','_pyfai.edf')
        mask_avGain = np.where(avimGain <0, 1, 0)
        poni.integrate1d(data = avimGain, filename = outfileGain,mask = mask_avGain,polarization_factor = 0.99,unit = '2th_deg',
                        correctSolidAngle = True, method = 'bbox',npt = 5000, error_model = 'poisson', safe = False)
        
        pattern2d, tth, eta, e = poni.integrate2d(data = avimGain, filename = outfileGain2d,mask = mask_avGain,polarization_factor = 0.99,unit = '2th_deg',
                        correctSolidAngle = True, method = 'bbox',npt_rad = 5000, npt_azim = 360, error_model = 'poisson', safe = False)
        bubbleHeader(outfileGain2d, pattern2d, tth, eta,)
        clearPyFAI_header(outfile)
        clearPyFAI_header(outfileGain)
        
        for c,file in enumerate(files):
            print(file)
            xyefile = file.replace('.cbf','.xye')
            outputfile = f'{dest}/{subdir}/{xyefile}'
            poni.integrate1d(data = dataset[:,:,c], filename = outputfile,mask = maskdct[c],polarization_factor = 0.99,unit = '2th_deg',
                            correctSolidAngle = True, method = 'bbox',npt = 5000, error_model = 'poisson', safe = False)
            dataGainCorr = gainCorrection(dataset[:,:,c],gainArray)
            maskGain = np.where(dataGainCorr < 0, 1, maskdct[c])
            outputfileGain = f'{dest}/{subdirGain}/{xyefile}'
            poni.integrate1d(data = dataGainCorr, filename = outputfileGain,mask = maskGain,polarization_factor = 0.99,unit = '2th_deg',
                            correctSolidAngle = True, method = 'bbox',npt = 5000, error_model = 'poisson', safe = False)
            clearPyFAI_header(outputfileGain)
            x,y,e = np.loadtxt(outputfile,unpack = True,comments = '#')
            xg,yg,eg = np.loadtxt(outputfileGain,comments = '#',unpack = True)
            if clearHeader:
                np.savetxt(outputfile,np.array([x,y,e]).transpose(),fmt = '%.6e')
            if c == 0:
                av1d = np.empty(shape = (len(y),len(files)))
                eav = np.empty(shape = (len(y),len(files)))
                av1dg = np.empty(shape = (len(yg),len(files)))
                eavg = np.empty(shape = (len(yg),len(files)))
            av1d[:,c] = y
            eav[:,c] = e
            av1dg[:,c] = yg
            eavg[:,c] = eg
        av1d = np.average(av1d,axis=1)
        av1dg = np.average(av1dg,axis=1)
        eav = np.average(eav,axis = 1)
        eavg = np.average(eavg,axis = 1)
        np.savetxt(f'{dest}/{subdir}/{avdir}/{xyefile}',np.array([x,av1d,eav]).transpose())
        np.savetxt(f'{dest}/{subdirGain}/{avdir}/{xyefile}',np.array([xg,av1dg,eavg]).transpose())
        
        qav = av1d*4*np.pi*np.sin(x*np.pi/(180*2))/wavelength
        np.savetxt(f'{dest}/{subdir}/{avdir}/qav_{i}.xy',np.array([x,qav]).transpose())
if __name__ == '__main__':
    run(direc,poni,mask,gainFile,averaging)