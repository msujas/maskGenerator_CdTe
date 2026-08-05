import fabio
import fabio.edfimage
import numpy as np
import matplotlib.pyplot as plt
import os, re
from cryio.cbfimage import CbfImage, CbfHeader
from glob import glob
import pyFAI
import argparse
from typing import Literal
import math

def gainCorrection(avim,gainArray):
    avimGain = avim/gainArray
    avimGain = np.where(np.bitwise_or(gainArray <0, avim<0), -1, avimGain)
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
    
def bubbleHeader(file2d,array2d, tth, eta, y, e):
    xye = np.array([tth,y,e]).transpose().flatten()
    xyestring = ' '.join([str(i) for i in xye])
    header = {
    'Bubble_cake_version' : 3,
    'Bubble_cake' : f'{tth[0]} {tth[-1]} {eta[0]} {eta[-1]}',
    'Bubble_normalized': 1 ,
    'Bubble_pattern': xyestring
    }
    f = fabio.edfimage.EdfImage(data = array2d[::-1,:], header = header)
    f.write(file2d)

def getAvFiles(direc, outdir = 'average'):
    files= glob(f'{direc}/*.cbf')
    basefilename = os.path.basename(files[-1])
    shortbasename = re.sub('_[0-9][0-9][0-9][0-9]p','',basefilename).replace('.cbf','')
    return f'{direc}/{outdir}/{shortbasename}_average.cbf', f'{direc}/{outdir}/{shortbasename}_average_gainCorrected.cbf'

class Integrator():
    def __init__(self, direc, maskfile, ponifile, gainfile=None, polfactor=0.85, fileextension='cbf', unit:Literal['2th_deg','q_A']='2th_deg',dest=None,
                 nstdevs = 3, avdir = 'average'):
        self.direc = direc
        self.badframeslog = f'{direc}/badframes.log'
        self.files = glob(f'{direc}/*.{fileextension}')
        self.basemask = fabio.open(maskfile).data
        if gainfile:
            self.flatfield = fabio.open(gainfile).data
        else:
            self.flatfield = None
        self.polfactor = polfactor
        self.unit = unit
        if not dest:
            dest = direc
        self.dest = dest
        self.nstdevs = nstdevs
        self.avdir = avdir
        self.ponifile = ponifile
        self.fileextension = fileextension
        

        self.ffcheck = not self.flatfield is None

    def appendBadFrames(self,file):
            f = open(self.badframeslog,'a')
            f.write(f'{file}\n')
            f.close()

    def makeDataSet(self, scale = 10**9, doMonitor = True):
        self.usedFiles = []
        monitors = []
        i1:np.ndarray = fabio.open(self.files[0]).data
        count = 0
        for file in self.files:
            header = CbfHeader(file)
            try:
                monitor = header['Flux']
            except:
                print(f'{file} flux not recorded, not including in averaging')
                self.appendBadFrames(self.badframeslog,file)
                continue
            exposure = header['Exposure_time']
            if monitor < 1000 * exposure:
                print(f'{file} low flux, not including in averaging')
                self.appendBadFrames(self.badframeslog,file)
                continue
            self.usedFiles.append(file)
            monitors.append(monitor)

        self.dataset = np.empty(shape = (*i1.shape,len(self.usedFiles)))

        for count,(file,monitor) in enumerate(zip(self.usedFiles,monitors)):
            cbf = CbfImage(file)
            array = cbf.array
            if doMonitor:
                array = (array/monitor)*scale
            else:
                array = array*1000 #multiply by 1000 as 10^6 is common monitor value
            self.dataset[:,:,count] = array

    def makeMasks(self,  plot = False):
        self.maskdct = {}
        #maskarray = np.zeros(shape=dataset.shape)
        median = np.median(self.dataset,axis=2)
        stdev = np.std(self.dataset,axis = 2)
        for c,file in enumerate(self.usedFiles):
            print(file)
            array = self.dataset[:,:,c]
            self.maskdct[c] = np.where(array > median+self.nstdevs*stdev,1,self.basemask)
            #maskarray[:,:,c] = np.where(array > median+nstdevs*stdev,1,baseMask)
            if plot:
                vmax = np.percentile(np.where(np.isnan(array),0,array),99.5)
                fig,ax = plt.subplots(1,2,dpi = 150)
                ax[0].imshow(self.maskdct[c])
                ax[1].imshow(array,vmax = vmax)
                plt.show()



    def integrateAverage(self,  npt = 5000, nptA = 360,  shortbasename = None):
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
        self.poni = pyFAI.load(self.ponifile)
        if not shortbasename:
            basefilename = os.path.basename(self.usedFiles[-1])
            shortbasename = re.sub('_[0-9][0-9][0-9][0-9]p','',basefilename).replace('.cbf','')

        if not os.path.exists(f'{self.dest}/{self.avdir}/xye/'):
            os.makedirs(f'{self.dest}/{self.avdir}/xye/')

        dataset2 = np.empty(shape = self.dataset.shape)
        for n in self.maskdct:
            dataset2[:,:,n] = np.where(self.maskdct[n] == 0, self.dataset[:,:,n], np.nan)
        avim = np.nanmean(dataset2, axis = 2)
        avim = np.where(np.isnan(avim), -2, avim)
        im = CbfImage()
        im.array = avim
        im.save(f'{self.dest}/{self.avdir}/{shortbasename}_average.cbf')
        
        if self.ffcheck:
            avimGain = im.array = gainCorrection(avim,self.flatfield)
            im.save(f'{self.dest}/{self.avdir}/{shortbasename}_average_gainCorrected.cbf')
            mask_avGain = np.where(avimGain < 0, 1, 0)


        outfile = f'{self.dest}/{self.avdir}/xye/{shortbasename}_average.xye'
        outfile_2d = outfile.replace('.xye','.edf')
        mask_av = np.where(avim < 0, 1, 0)
        
        x,y,e = self.poni.integrate1d(data = avim, filename = outfile,mask = mask_av,polarization_factor = self.polfactor,unit = self.unit,
                        correctSolidAngle = True, method = 'bbox',npt = npt, error_model = 'poisson', safe = False)
        print(f'saving {outfile}')
        np.savetxt(outfile,np.array([x,y,e]).transpose(), fmt="%.6f")
        print(f'saving {outfile_2d}')
        result = self.poni.integrate2d(data = avim, filename = outfile_2d,mask = mask_av,polarization_factor = self.polfactor,unit = self.unit,
                        correctSolidAngle = True, method = 'bbox',npt_rad = npt, npt_azim = nptA, error_model = 'poisson', safe = False)
        bubbleHeader(outfile_2d,*result[:3], y, e)

        if self.ffcheck:
            outfileGC = f'{self.dest}/{self.avdir}/xye/{shortbasename}_average_gainCorrected.xye'
            print(f'saving {outfileGC}')
            x,y,e = self.poni.integrate1d(data = avimGain, filename = outfileGC,mask =mask_avGain,polarization_factor = self.polfactor,unit = self.unit,
                            correctSolidAngle = True, method = 'bbox',npt = npt, error_model = 'poisson', safe = False)
            outfileGC_2d = outfileGC.replace('.xye','.edf')
            print(f'saving {outfileGC_2d}')
            result = self.poni.integrate2d(data = avimGain, filename = outfileGC_2d,mask = mask_avGain,polarization_factor = self.polfactor,unit = self.unit,
                            correctSolidAngle = True, method = 'bbox',npt_rad = npt,npt_azim = nptA, error_model = 'poisson', safe = False)
            clearPyFAI_header(outfileGC)
            bubbleHeader(outfileGC_2d,*result[:3], y, e)

    def integrateIndividual(self, subdir,  npt = 5000):
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
        self.poni = pyFAI.load(self.ponifile)
        while subdir[-1] == '/' or subdir[-1] == '\\':
            subdir = subdir[:-1]
        subdirGain = f'{subdir}_gainCorr'
        wavelength = self.poni.wavelength*10**10
        if not os.path.exists(f'{self.dest}/{subdir}/{self.avdir}/'):
            os.makedirs(f'{self.dest}/{subdir}/{self.avdir}/')
        if not os.path.exists(f'{self.dest}/{subdirGain}/{self.avdir}/') and self.flatfield:
            os.makedirs(f'{self.dest}/{subdirGain}/{self.avdir}/')  
        for c,file in enumerate(self.usedFiles):
            print(file)
            xyefile = file.replace('.cbf','.xye')
            outputfile = f'{self.dest}/{subdir}/{xyefile}'
            print(f'saving {outputfile}')
            x,y,e = self.poni.integrate1d(data = self.dataset[:,:,c], filename = outputfile,mask = self.maskdct[c],polarization_factor = self.polfactor,unit = self.unit,
                            correctSolidAngle = True, method = 'bbox',npt = npt, error_model = 'poisson', safe = False)
            np.savetxt(outputfile, np.array([x,y,e]).transpose(),fmt = "%.6f")
            if self.ffcheck:
                outputfileGC = f'{self.dest}/{subdirGain}/{xyefile}'
                print(f'saving {outputfileGC}')
                arrayGC = gainCorrection(self.dataset[:,:,c],self.flatfield)
                xg,yg,eg = self.poni.integrate1d(data = arrayGC, filename = outputfileGC,mask = self.maskdct[c],polarization_factor = self.polfactor,unit = self.unit,
                                correctSolidAngle = True, method = 'bbox',npt = npt, error_model = 'poisson', safe = False)
                clearPyFAI_header(outputfileGC)

            if c == 0:
                av1d = np.empty(shape = (len(y),len(self.usedFiles)))
                eav = np.empty(shape = (len(y),len(self.usedFiles)))
                if self.ffcheck:
                    av1dg = np.empty(shape = (len(yg),len(self.usedFiles)))
                    eavg = np.empty(shape = (len(yg),len(self.usedFiles)))
            av1d[:,c] = y
            eav[:,c] = e
            if self.ffcheck:
                av1dg[:,c] = yg
                eavg[:,c] = eg
        av1d = np.average(av1d,axis=1)
        eav = np.average(eav,axis = 1)
        if self.ffcheck:
            av1dg = np.average(av1dg,axis=1)
            eavg = np.average(eavg,axis = 1)
        np.savetxt(f'{self.dest}/{subdir}/average/{xyefile}',np.array([x,av1d,eav]).transpose())

        qav = av1d*4*np.pi*np.sin(x*np.pi/(180*2))/wavelength
        np.savetxt(f'{self.dest}/{subdir}/average/qav.xy',np.array([x,qav]).transpose())
        
        if self.ffcheck:
            qavg = av1dg*4*np.pi*np.sin(xg*np.pi/(180*2))/wavelength
            np.savetxt(f'{self.dest}/{subdirGain}/average/{xyefile}',np.array([xg,av1dg,eavg]).transpose())
            np.savetxt(f'{self.dest}/{subdirGain}/average/qav.xy',np.array([xg,qavg]).transpose())
        
    
    def runAll(self, scale = 10**9, individual=False, plotmasks = False):
        self.makeDataSet(scale=scale)
        self.makeMasks(plotmasks)
        self.integrateAverage()
        if individual:
            self.integrateIndividual()

    def runsplit(self, split, scale=1e9,individual=False, plotmasks=False):
        allfiles = glob(f'{self.direc}/*.{self.fileextension}')
        if not split:
            split = len(allfiles)
        for i in range(math.ceil(len(allfiles)/split)):
            files = allfiles[i*split:i*split+split]
            self.files = files
            self.runAll(scale,individual,plotmasks)


def averagefiles(*files, basemask, flatfieldfile=None):
    filelist = []
    for f in files:
        filelist += glob(f)
    dirname = os.path.dirname(os.path.abspath(filelist[0]))
    i = Integrator('',maskfile=basemask, gainfile=flatfieldfile, ponifile='')
    i.files = filelist
    i.makeDataSet()
    i.makeMasks()
    for n in range(i.dataset.shape[2]):
        i.dataset[:,:,n] = np.where(i.maskdct[n] >0, np.nan, i.dataset[:,:,n])
    av = np.nanmean(i.dataset, axis=2)
    av = np.where(np.isnan(av),-1, av)
    avdir = f'{dirname}/average'
    os.makedirs(avdir,exist_ok=True)
    im = CbfImage()
    im.array = av
    im.header['Flux'] = 1
    basefilename = os.path.basename(filelist[-1])
    shortbasename = re.sub('_[0-9][0-9][0-9][0-9]p','',basefilename).replace('.cbf','')
    outfile = f'{avdir}/{shortbasename}.cbf'
    print(f'saving {outfile}')
    im.save_cbf(outfile)
    if flatfieldfile:
        ff:np.ndarray = fabio.open(flatfieldfile).data
        avgc = gainCorrection(av,ff)
        imgc = CbfImage()
        imgc.array = avgc
        imgc.header['Flux'] = 1
        outfileGC = f'{avdir}/{shortbasename}GC.cbf'
        print(f'saving {outfileGC}')
        imgc.save(outfileGC)
    return av

def averagefilescli():
    ap = argparse.ArgumentParser()
    ap.add_argument('-m','--mask', type=str, default=None, help='maskfile used for averageing')
    ap.add_argument('files', nargs='+', help = 'files to average')
    ap.add_argument('-ff','--flatfield', type=str, default=None, help='flat field (gain map) file used to correct image')
    args = ap.parse_args()
    mask = args.mask
    files = args.files
    fffile = args.flatfield
    averagefiles(*files,basemask= mask,flatfieldfile= fffile)