import fabio
import fabio.edfimage
import numpy as np
import matplotlib.pyplot as plt
import os, re
from cryio.cbfimage import CbfImage, CbfHeader

def gainCorrection(avim,gainArray):
    avimGain = avim/gainArray
    avimGain = np.where(gainArray <0, -1, avimGain)
    return avimGain

def gainCorrectionFiles(cbfFile, gainFile):
    array = CbfImage(cbfFile).array
    gainArray = fabio.open(gainFile).data
    newarray = gainCorrection(array,gainArray)
    im = fabio.edfimage.EdfImage(newarray)
    im.save(cbfFile.replace('.cbf','_GC.edf'))

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

def appendBadFrames(badFramesLog,file):
        f = open(badFramesLog,'a')
        f.write(f'{file}\n')
        f.close()

def makeDataSet(files : list, badFramesLog : str, scale = 10**9, doMonitor = True):
    usedFiles = []
    monitors = []
    i1 = fabio.open(files[0]).data
    count = 0
    for file in files:
        header = CbfHeader(file)
        try:
            monitor = header['Flux']
        except:
            print(f'{file} flux not recorded, not including in averaging')
            appendBadFrames(badFramesLog,file)
            continue
        exposure = header['Exposure_time']
        if monitor < 1000 * exposure:
            print(f'{file} low flux, not including in averaging')
            appendBadFrames(badFramesLog,file)
            continue
        usedFiles.append(file)
        monitors.append(monitor)

    dataset = np.empty(shape = (*i1.shape,len(usedFiles)))

    for count,(file,monitor) in enumerate(zip(usedFiles,monitors)):
        cbf = CbfImage(file)
        array = cbf.array
        if doMonitor:
            array = (array/monitor)*scale
        else:
            array = array*1000 #multiply by 1000 as 10^6 is common monitor value
        dataset[:,:,count] = array
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

def integrateAverage(dataset, files, dest, poni, gainFile, maskdct, outdir='average', unit = '2th_deg', npt = 5000, nptA = 360, polF = 0.99, shortbasename = None):
    '''
    arguments: dataset, files, dest, poni, gainFile, maskdct, unit = '2th_deg', npt = 5000, nptA = 360, polF = 0.99
    dataset - array containing individual diffraction images, shape = (y image size, x image size, number of images)
    files - files used to make dataset
    dest - destination directory
    poni - pyfai azimuthal integrator object
    gainFile - file where gain array is stored
    maskdct - dictionary of masks, same length as dataset
    unit - integration unit ('2th_deg', 'q_A', 'q_nm')
    npt - number of radial points
    nptA - number of azimuthal points
    polF - polarisation factor

    creates and average image, with each image masked individually. Integrates the average with and without gain correction
    no returned variable
    '''
    if shortbasename == None:
        basefilename = os.path.basename(files[-1])
        shortbasename = re.sub('_[0-9][0-9][0-9][0-9]p','',basefilename).replace('.cbf','')

    if not os.path.exists(f'{dest}/{outdir}/xye/'):
        os.makedirs(f'{dest}/{outdir}/xye/')

    dataset2 = np.empty(shape = dataset.shape)
    for n in maskdct:
        dataset2[:,:,n] = np.where(maskdct[n] == 0, dataset[:,:,n], np.nan)
    avim = np.nanmean(dataset2, axis = 2)
    avim = np.where(np.isnan(avim), -2, avim)
    im = CbfImage()
    im.array = avim
    im.save(f'{dest}/{outdir}/{shortbasename}_average.cbf')
    
    if gainFile != None:
        gainArray = fabio.open(gainFile).data
        avimGain = im.array = gainCorrection(avim,gainArray)
        im.save(f'{dest}/{outdir}/{shortbasename}_average_gainCorrected.cbf')
        mask_avGain = np.where(avimGain < 0, 1, 0)


    outfile = f'{dest}/{outdir}/xye/{shortbasename}_average.xye'
    outfile_2d = outfile.replace('.xye','.edf')
    mask_av = np.where(avim < 0, 1, 0)
    
    x,y,e = poni.integrate1d(data = avim, filename = outfile,mask = mask_av,polarization_factor = polF,unit = unit,
                    correctSolidAngle = True, method = 'bbox',npt = npt, error_model = 'poisson', safe = False)
    np.savetxt(outfile,np.array([x,y,e]).transpose(), fmt="%.6f")
    result = poni.integrate2d(data = avim, filename = outfile_2d,mask = mask_av,polarization_factor = polF,unit = unit,
                    correctSolidAngle = True, method = 'bbox',npt_rad = npt, npt_azim = nptA, error_model = 'poisson', safe = False)
    bubbleHeader(outfile_2d,*result[:3])

    if gainFile != None:
        outfileGC = f'{dest}/average/xye/{shortbasename}_average_gainCorrected.xye'
        poni.integrate1d(data = avimGain, filename = outfileGC,mask =mask_avGain,polarization_factor = polF,unit = unit,
                        correctSolidAngle = True, method = 'bbox',npt = npt, error_model = 'poisson', safe = False)
        outfileGC_2d = outfileGC.replace('.xye','.edf')
        result = poni.integrate2d(data = avimGain, filename = outfileGC_2d,mask = mask_avGain,polarization_factor = polF,unit = unit,
                        correctSolidAngle = True, method = 'bbox',npt_rad = npt,npt_azim = nptA, error_model = 'poisson', safe = False)
        clearPyFAI_header(outfileGC)
        bubbleHeader(outfileGC_2d,*result[:3])

def integrateIndividual(dataset,files, dest, subdir, poni, maskdct, gainFile, avdir = 'average', unit = '2th_deg', npt = 5000, polF = 0.99):
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
    if not os.path.exists(f'{dest}/{subdirGain}/{avdir}/') and gainFile != None:
        os.makedirs(f'{dest}/{subdirGain}/{avdir}/')  
    for c,file in enumerate(files):
        print(file)
        xyefile = file.replace('.cbf','.xye')
        outputfile = f'{dest}/{subdir}/{xyefile}'
        x,y,e = poni.integrate1d(data = dataset[:,:,c], filename = outputfile,mask = maskdct[c],polarization_factor = polF,unit = unit,
                        correctSolidAngle = True, method = 'bbox',npt = npt, error_model = 'poisson', safe = False)
        np.savetxt(outputfile, np.array([x,y,e]).transpose(),fmt = "%.6f")
        if gainFile != None:
            gainArray = fabio.open(gainFile).data
            outputfileGC = f'{dest}/{subdirGain}/{xyefile}'
            arrayGC = dataset[:,:,c]/gainArray
            arrayGC = np.where(gainArray < 0, -1, arrayGC)
            xg,yg,eg = poni.integrate1d(data = arrayGC, filename = outputfileGC,mask = maskdct[c],polarization_factor = polF,unit = unit,
                            correctSolidAngle = True, method = 'bbox',npt = npt, error_model = 'poisson', safe = False)
            clearPyFAI_header(outputfileGC)

        
        


        if c == 0:
            av1d = np.empty(shape = (len(y),len(files)))
            eav = np.empty(shape = (len(y),len(files)))
            if gainFile != None:
                av1dg = np.empty(shape = (len(yg),len(files)))
                eavg = np.empty(shape = (len(yg),len(files)))
        av1d[:,c] = y
        eav[:,c] = e
        if gainFile != None:
            av1dg[:,c] = yg
            eavg[:,c] = eg
    av1d = np.average(av1d,axis=1)
    eav = np.average(eav,axis = 1)
    if gainFile != None:
        av1dg = np.average(av1dg,axis=1)
        eavg = np.average(eavg,axis = 1)
    np.savetxt(f'{dest}/{subdir}/average/{xyefile}',np.array([x,av1d,eav]).transpose())

    qav = av1d*4*np.pi*np.sin(x*np.pi/(180*2))/wavelength
    np.savetxt(f'{dest}/{subdir}/average/qav.xy',np.array([x,qav]).transpose())
    
    if gainFile != None:
        qavg = av1dg*4*np.pi*np.sin(xg*np.pi/(180*2))/wavelength
        np.savetxt(f'{dest}/{subdirGain}/average/{xyefile}',np.array([xg,av1dg,eavg]).transpose())
        np.savetxt(f'{dest}/{subdirGain}/average/qav.xy',np.array([xg,qavg]).transpose())
    
    