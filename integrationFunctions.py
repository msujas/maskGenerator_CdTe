import fabio
import numpy as np
import matplotlib.pyplot as plt
import os, re
from cryio.cbfimage import CbfImage

def gainCorrection(avim,gainArray):
    avimGain = avim/gainArray
    avimGain = np.where(gainArray <0, -1, avimGain)
    return avimGain

def clearPyFAI_header(file):
    x,y,e = np.loadtxt(file,unpack = True, comments = '#')
    np.savetxt(file,np.array([x,y,e]).transpose(), '%.6f')
    
def bubbleHeader(file2d,array2d, tth, eta):
    header = {
    'Bubble_cake' : f'{tth[0]} {tth[-1]} {eta[0]} {eta[-1]}',
    'Bubble_normalized': 1 
    }
    f = fabio.edfimage.EdfImage(data = array2d.transpose(), header = header)

    f.write(file2d)

def makeDataSet(files : list, badFramesLog : str, scale = 10**9, doMonitor = True):
    usedFiles = []
    i1 = fabio.open(files[0]).data
    dataset = np.empty(shape = (*i1.shape,len(files)))
    count = 0
    for file in files:
        cbf = CbfImage(file)
        array = cbf.array
        header = cbf.header
         
        if doMonitor:
            try:
                monitor = header['Flux']
                exposure = header['Exposure_time']

                if monitor <= 1000*exposure:
                    f = open(badFramesLog,'a')
                    f.write(f'{file}\n')
                    f.close()
                    print(f'{file} low flux, not including in averaging')
                    dataset = dataset[:,:,:-1]
                    continue
            except:
                f = open(badFramesLog,'a')
                f.write(f'{file}\n')
                f.close()
                print(f'{file} flux not recorded, not including in averaging')
                dataset = dataset[:,:,:-1]
                continue
            array = (array/monitor)*scale
        else:
            array = array*1000 #multiply by 1000 as 10^6 is common monitor value
        dataset[:,:,count] = array
        usedFiles.append(file)
        count += 1
            
    return dataset, usedFiles

def makeMasks(dataset, files, baseMask, nstdevs = 3, plot = False):
    maskdct = {}
    median = np.median(dataset,axis=2)
    stdev = np.std(dataset,axis = 2)
    for c,file in enumerate(files):
        print(file)
        array = dataset[:,:,c]
        
        maskdct[c] = np.where(array > median+nstdevs*stdev,1,baseMask)
        if plot:
            vmax = np.percentile(np.where(np.isnan(array),0,array),99.5)
            fig,ax = plt.subplots(1,2,dpi = 150)
            ax[0].imshow(maskdct[c])
            ax[1].imshow(array,vmax = vmax)
            plt.show()
    return maskdct

def integrateAverage(dataset, files, dest, poni, gainArray, maskdct, unit = '2th_deg', npt = 5000, nptA = 360, polF = 0.99, fileAppend = ''):
    '''
    arguments: dataset, files, dest, poni, gainArray, maskdct, unit = '2th_deg', npt = 5000, nptA = 360, polF = 0.99
    dataset - array containing individual diffraction images, shape = (y image size, x image size, number of images)
    files - files used to make dataset
    dest - destination directory
    poni - pyfai azimuthal integrator object
    gainArray - array containing gain values for images
    maskdct - dictionary of masks, same length as dataset
    unit - integration unit ('2th_deg', 'q_A', 'q_nm')
    npt - number of radial points
    nptA - number of azimuthal points
    polF - polarisation factor

    creates and average image, with each image masked individually. Integrates the average with and without gain correction
    no returned variable
    '''
    basefilename = os.path.basename(files[-1])
    shortbasename = re.sub('_[0-9][0-9][0-9][0-9]p',fileAppend,basefilename).replace('.cbf','')

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
    poni.integrate1d(data = avim, filename = outfile,mask = mask_av,polarization_factor = polF,unit = unit,
                    correctSolidAngle = True, method = 'bbox',npt = npt, error_model = 'poisson', safe = False)
    result = poni.integrate2d(data = avim, filename = outfile_2d,mask = mask_av,polarization_factor = polF,unit = unit,
                    correctSolidAngle = True, method = 'bbox',npt_rad = npt, npt_azim = nptA, error_model = 'poisson', safe = False)
    bubbleHeader(outfile_2d,*result[:3])

    outfileGC = f'{dest}/average/xye/{shortbasename}_average_gainCorrected.xye'
    poni.integrate1d(data = avimGain, filename = outfileGC,mask =mask_avGain,polarization_factor = polF,unit = unit,
                    correctSolidAngle = True, method = 'bbox',npt = npt, error_model = 'poisson', safe = False)
    outfileGC_2d = outfileGC.replace('.xye','_pyfai.edf')
    result = poni.integrate2d(data = avimGain, filename = outfileGC_2d,mask = mask_avGain,polarization_factor = polF,unit = unit,
                    correctSolidAngle = True, method = 'bbox',npt_rad = npt,npt_azim = nptA, error_model = 'poisson', safe = False)
    clearPyFAI_header(outfile)
    clearPyFAI_header(outfileGC)
    bubbleHeader(outfileGC_2d,*result[:3])


def integrateIndividual(dataset,files, dest, subdir, poni, maskdct, gainArray, avdir = 'average', unit = '2th_deg', npt = 5000, polF = 0.99):
    '''
    argumenets: dataset,files, dest, subdir, poni, maskdct, gainArray, avdir = 'average', unit = '2th_deg', npt = 5000, polF = 0.99
    dataset - array containing individual diffraction images, shape = (y image size, x image size, number of images)
    files - files used to make dataset
    dest - destination directory
    subdir - subdirectory where individual xye files are saved
    poni - pyfai azimuthal integrator object
    maskdct - dictionary of masks, same length as dataset
    gainArray - array containing gain values for images
    avdir - name of directory where the average of the integrated files is saved
    unit - integration unit ('2th_deg', 'q_A', 'q_nm')
    npt - number of radial points
    polF - polarisation factor

    integrates each image in dataset with it's own mask, with and without gain correction, then averages the integrated patterns at the end.
    '''
    while subdir[-1] == '/' or subdir[-1] == '\\':
        subdir = subdir[:-1]
    subdirGain = f'{subdir}_gainCorr'
    wavelength = poni.wavelength*10**10
    if not os.path.exists(f'{dest}/{subdir}/{avdir}/'):
        os.makedirs(f'{dest}/{subdir}/{avdir}/')
    if not os.path.exists(f'{dest}/{subdirGain}/{avdir}/'):
        os.makedirs(f'{dest}/{subdirGain}/{avdir}/')  
    for c,file in enumerate(files):
        print(file)
        xyefile = file.replace('.cbf','.xye')
        outputfile = f'{dest}/{subdir}/{xyefile}'
        x,y,e = poni.integrate1d(data = dataset[:,:,c], filename = outputfile,mask = maskdct[c],polarization_factor = polF,unit = unit,
                        correctSolidAngle = True, method = 'bbox',npt = npt, error_model = 'poisson', safe = False)

        outputfileGC = f'{dest}/{subdirGain}/{xyefile}'
        arrayGC = dataset[:,:,c]/gainArray
        arrayGC = np.where(gainArray < 0, -1, arrayGC)
        xg,yg,eg = poni.integrate1d(data = arrayGC, filename = outputfileGC,mask = maskdct[c],polarization_factor = polF,unit = unit,
                        correctSolidAngle = True, method = 'bbox',npt = npt, error_model = 'poisson', safe = False)
        

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