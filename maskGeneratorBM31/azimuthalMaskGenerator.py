import fabio.edfimage
import numpy as np
import fabio
from cryio.cbfimage import CbfHeader, CbfImage
import matplotlib.pyplot as plt
import pyFAI
import pyFAI.geometry
import os
from glob import glob
if __name__ == '__main__':
   from integrationFunctions import clearPyFAI_header, gainCorrection, bubbleHeader
else:
    from . import clearPyFAI_header, gainCorrection


datadir = r''
poni = r''
maskfile = r''
polarisation = 0.99
stdevs = 4
threshold = 100 #counts above 1 stdev to mask
scale = 10**5
gainFile = None #r'C:\Users\kenneth1a\Documents\beamlineData\gainMaps/gainMap_filtered_kpm_2025-02-05.edf'
saveMasks = True
save2d = True
nbins = 800

def generateMask(dataarray, binarray,nbins, stdevs, basemask, threshold = 100):
    dataarray2 = np.where(dataarray < 0, np.nan, dataarray)
    array = np.empty(shape = dataarray.shape)
    newmask = basemask.copy()
    for n in range(nbins):
        wherebin = np.where(binarray == n)
        array =  dataarray2[wherebin]
        maskBin = basemask[wherebin]
        stdev = np.nanstd(array)
        median = np.nanmedian(array)
        maskBin = np.where((array > median + stdevs*stdev) | (array > median + stdev + threshold) , 1, maskBin)
        y = np.where(binarray==n)[0]
        x = np.where(binarray==n)[1]
        newmask[[y],[x]] = maskBin
    return newmask

def makeArray(file, polArray, saArray, gainMap = None):
    dataarray = CbfImage(file).array
    dataarray = dataarray/polArray
    dataarray = dataarray/ saArray
    if type(gainMap) == np.ndarray:
        dataarray = gainCorrection(dataarray, gainMap)
    return dataarray

def generateDetArrays(cbffile, ponifile, nbins):
    dataarray = CbfImage(cbffile).array
    poni = pyFAI.load(ponifile)
    geometry = pyFAI.geometry.Geometry()
    geometry.load(poni)
    array2th = geometry.twoThetaArray()*180/np.pi
    polarray  = geometry.polarization(factor = polarisation)
    saArray = geometry.solidAngleArray()
    binsize = np.max(array2th)/nbins
    binarray = (array2th/binsize).astype(np.uint16)
    binarray = np.where(dataarray < 0, -1, binarray)
    return array2th, polarray, saArray, binarray

def integrateIndividualAzMask(file, maskfile, ponifile, gainFile=None, stdevs = 3, nbins = 800, outdir = 'xye'):
    direc = os.path.dirname(file)
    fulloutdir = f'{direc}/{outdir}/'
    outbasefile = os.path.basename(file).replace('.cbf','.xye')
    outfile = f'{fulloutdir}/{outbasefile}'
    array2th, polarray, saArray, binarray = generateDetArrays(file, ponifile, nbins)
    if not gainFile == None:
        gainMap = fabio.open(gainMap).data
    else:
        gainMap = None
    basemask = fabio.open(maskfile).data
    dataarray = makeArray(file,polarray,saArray, gainMap)
    mask = generateMask(dataarray, binarray, nbins, stdevs, basemask)
    poni = pyFAI.load(ponifile)
    poni.integrate1d(dataarray, mask = mask, filename=outfile, polarization_factor=0.99, unit = '2th_deg', correctSolidAngle=True,
                     method = 'bbox', npt = 5000, error_model='poisson', safe = False)

def run(datadir, ponifile,  stdevs, maskfile, scale, threshold = 100, polarisation = 0.99, gainFile = None, nbins = 800, ext = 'cbf',
        outdir = 'xye', save2d = False, saveMasks = False):
    os.chdir(datadir)
    if not os.path.exists(f'{datadir}/{outdir}/'):
        os.makedirs(f'{datadir}/{outdir}/')
    basemask = fabio.open(maskfile).data
    cbfs = glob(f'*.{ext}')
    array2th, polarray, saArray, binarray = generateDetArrays(cbfs[0],ponifile,nbins)
    poni = pyFAI.load(ponifile)
    maskdir = f'{datadir}/masks{stdevs}/'
    if gainFile != None:
        gainMap = fabio.open(gainFile).data
    else:
        gainMap=None
    os.makedirs(maskdir,exist_ok=True)
    for file in cbfs:
        print(file)
        dataarray = makeArray(file,polarray, saArray,gainMap)
        xyefile = file.replace('.cbf','.xye')
        outfile = f'{outdir}/{xyefile}'
        header = CbfHeader(file)
        monitorCounts = header['Flux']
        mask = generateMask(dataarray,binarray,nbins, stdevs, basemask, threshold)
        if saveMasks:
            maskim = fabio.edfimage.EdfImage(mask)
            maskim.save(f'{maskdir}/{file}'.replace('.cbf','.edf'))
        normArray = (dataarray/monitorCounts) * scale
        poni.integrate1d(normArray,mask = mask, filename = outfile, polarization_factor = polarisation,
                        unit = '2th_deg', correctSolidAngle = True, method = 'bbox', npt = 5000, 
                        error_model = 'poisson', safe = False)
        clearPyFAI_header(outfile)
        if save2d:
            outfile_2d = outfile.replace('.xye','.edf')
            result = poni.integrate2d(data = normArray, filename = outfile_2d,mask = mask,polarization_factor = polarisation, unit = '2th_deg',
                        correctSolidAngle = True, method = 'bbox',npt_rad = 5000,npt_azim = 360, error_model = 'poisson', safe = False)
            bubbleHeader(outfile_2d,*result[:3])
            
def runRecursive(direc, ponifile, maskfile, polarisation = 0.99, gainFile=None, stdevs = 4, scale=1, threshold = 100, nbins= 800, 
                 ext = 'cbf', outdir = 'xye'):
    for root, dirs, files in os.walk(direc):
        cbfs = glob(f'{root}/*.{ext}')
        if not cbfs:
            continue
        run(root,ponifile, stdevs, maskfile, scale, threshold, polarisation, gainFile, nbins, ext,outdir, save2d=False, saveMasks=False)

if __name__ == '__main__':
    run(datadir, poni,  stdevs, maskfile, scale, threshold, polarisation,gainFile, save2d=save2d, saveMasks=saveMasks, nbins=nbins)