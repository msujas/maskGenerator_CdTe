import fabio.edfimage
import numpy as np
import fabio
from cryio.cbfimage import CbfHeader, CbfImage
import matplotlib.pyplot as plt
import pyFAI
import pyFAI.geometry
import os
from glob import glob
try:
   from integrationFunctions import clearPyFAI_header, gainCorrection, bubbleHeader
except ImportError:
    from . import clearPyFAI_header, gainCorrection, bubbleHeader
try:
    from azimuthCP import makeMaskCP
except ImportError:
    print('couldn\'t find cpp extension')

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

def generateMask(dataarray, basemask, binarray,nbins, stdevs,  threshold = 100, cpp = False):
    #print(dataarray.shape, basemask.shape, binarray.shape)
    if cpp:
        args = (dataarray.tolist(), basemask.tolist(), binarray.tolist(), nbins,stdevs,threshold)
        newmask = makeMaskCP(args)
        newmask = np.array(newmask)
        return newmask
    #python/numpy implementation
    dataarray2 = np.where(dataarray < 0, np.nan, dataarray)
    array = np.empty(shape = dataarray.shape)
    newmask = basemask.copy()
    for n in range(nbins):
        wherebin = np.where((binarray == n) & (basemask == 0))       
        if len(wherebin[0]) == 0:
            continue
        array =  dataarray2[wherebin]
        maskBin = basemask[wherebin]
        stdev = np.nanstd(array)
        median = np.nanmedian(array)
        maskBin = np.where((array > median + stdevs*stdev) | (array > median + stdev + threshold) , 1, maskBin)
        y = wherebin[0]
        x = wherebin[1]
        newmask[[y],[x]] = maskBin
    return newmask

def makeArray(file, polArray, saArray, gainMap = None):
    dataarray = CbfImage(file).array
    dataarray = dataarray/polArray
    dataarray = dataarray/ saArray
    if type(gainMap) == np.ndarray:
        dataarray = gainCorrection(dataarray, gainMap)
    return dataarray

def generateDetArrays(cbffile, ponifile, nbins, polarisation):
    dataarray = CbfImage(cbffile).array
    poni = pyFAI.load(ponifile)
    geometry = pyFAI.geometry.Geometry()
    geometry.load(poni)
    array2th = geometry.twoThetaArray()*180/np.pi
    polarray  = geometry.polarization(factor = polarisation)
    saArray = geometry.solidAngleArray()
    binsize = np.max(array2th)/nbins
    binarray = ((array2th+binsize/2)*(nbins-1)//np.max(array2th)).astype(np.uint16)
    binarray = np.where(dataarray < 0, -1, binarray)
    return array2th, polarray, saArray, binarray

def int2d(outfile1d, normArray, poni, mask):
    outfile_2d = outfile1d.replace('.xye','.edf')
    result = poni.integrate2d(data = normArray, filename = None,mask = mask,polarization_factor = None, unit = '2th_deg',
                correctSolidAngle = False, method = 'bbox',npt_rad = 5000,npt_azim = 360, error_model = 'poisson', safe = False)
    bubbleHeader(outfile_2d,*result[:3])
    return result

def integrateIndividualAzMask(file, maskfile, ponifile, gainFile=None, stdevs = 3, nbins = 800, outdir = 'xye', save2d = False):
    direc = os.path.dirname(file)
    fulloutdir = f'{direc}/{outdir}/'
    os.makedirs(fulloutdir,exist_ok=True)
    outbasefile = os.path.basename(file).replace('.cbf','.xye')
    outfile = f'{fulloutdir}/{outbasefile}'
    array2th, polarray, saArray, binarray = generateDetArrays(file, ponifile, nbins)
    if not gainFile == None:
        gainMap = fabio.open(gainFile).data
    else:
        gainMap = None
    basemask = fabio.open(maskfile).data.astype(np.uint8)
    dataarray = makeArray(file,polarray,saArray, gainMap)
    mask = generateMask(dataarray, basemask, binarray, nbins, stdevs)
    poni = pyFAI.load(ponifile)
    poni.integrate1d(dataarray, mask = mask, filename=outfile, polarization_factor=None, unit = '2th_deg', correctSolidAngle=False,
                     method = 'bbox', npt = 5000, error_model='poisson', safe = False)
    clearPyFAI_header(outfile)
    if save2d:
        int2d(outfile, dataarray, poni, mask)

def run(datadir, ponifile,  stdevs, maskfile, scale, threshold = 100, polarisation = 0.99, gainFile = None, nbins = 800, ext = 'cbf',
        outdir = 'xye', save2d = False, saveMasks = False, cpp = False, allFiles = None, averaging = 1):
    os.chdir(datadir)
    if not os.path.exists(f'{datadir}/{outdir}/'):
        os.makedirs(f'{datadir}/{outdir}/')
    outdirav = f'{datadir}/{outdir}_av{averaging}'
    if averaging > 1 and not os.path.exists(outdirav):
        os.makedirs(outdirav)
        
    basemask = fabio.open(maskfile).data
    cbfs = glob(f'*.{ext}')
    cbfs.sort()
    array2th, polarray, saArray, binarray = generateDetArrays(cbfs[0],ponifile,nbins, polarisation)
    poni = pyFAI.load(ponifile)
    maskdir = f'{datadir}/masks{stdevs}'
    if cpp:
        maskdir += 'cpp'
    if gainFile != None:
        gainMap = fabio.open(gainFile).data
    else:
        gainMap=None
    if saveMasks:
        os.makedirs(maskdir,exist_ok=True)
    if allFiles == None:
        allFiles = []
    allFilesR = [f for f in allFiles]
    av = 0

    for file in cbfs:
        fullfile = f'{datadir}/{file}'
        if fullfile not in allFilesR:
            allFilesR.append(fullfile)
        else:
            continue
        print(file)
        dataarray = makeArray(file,polarray, saArray,gainMap)
        xyefile = file.replace('.cbf','.xye')
        outfile = f'{outdir}/{xyefile}'
        header = CbfHeader(file)
        monitorCounts = 1
        if 'Flux' in header:
            monitorCounts = header['Flux']
        
        mask = generateMask(dataarray, basemask,binarray,nbins, stdevs, threshold, cpp)
        
        if saveMasks:
            maskim = fabio.edfimage.EdfImage(mask)
            maskim.save(f'{maskdir}/{file}'.replace('.cbf','.edf'))
        normArray = (dataarray/monitorCounts) * scale
        
        x,y,e = poni.integrate1d(normArray,mask = mask, filename = outfile, polarization_factor = None,
                        unit = '2th_deg', correctSolidAngle = False, method = 'bbox', npt = 5000, 
                        error_model = 'poisson', safe = False) #not applying polarisation and solid angle as already applied earlier
        if averaging > 1:
            if av == 0:
                yarray = np.empty(shape = (len(y),averaging))
                earray = np.empty(shape = (len(y),averaging))
            yarray[:,av] = y
            earray[:,av] = e
            if av == averaging - 1:
                outfileav = f'{outdirav}/{xyefile}'
                yarray = np.mean(yarray,axis = 1)
                earray = np.mean(earray,axis = 1)
                np.savetxt(outfileav, np.array([x,yarray,earray]).transpose(), fmt = "%.6f")
                
        clearPyFAI_header(outfile)
        print('file integrated')
        if save2d:
            result = int2d(outfile, normArray, poni, mask)
            print('cake saved')
            if averaging > 1:
                if av == 0:
                    resultav = np.empty(shape = (*result[0].shape, averaging))
                resultav[:,:,av] = result[0]
                if av == averaging - 1:
                    resultav = np.mean(resultav, axis = 2)
                    #cakeim = fabio.edfimage.EdfImage(resultav)
                    outfileav2d = outfileav.replace('.xye','.edf')
                    #cakeim.save(outfileav2d)
                    bubbleHeader(outfileav2d, resultav, result[1],result[2])
        if av == averaging -1:
            av = -1
        av += 1

            
    return allFilesR
            
def runRecursive(direc, ponifile, maskfile, polarisation = 0.99, gainFile=None, stdevs = 4, scale=1, threshold = 100, nbins= 800, 
                 ext = 'cbf', outdir = 'xye', cpp = False, saveMasks = False, save2d = False, allFiles = None, averaging : int = 1):
    if allFiles == None:
        allFiles = []
    allFilesR = [f for f in allFiles]
    for root, dirs, files in os.walk(direc):
        cbfs = glob(f'{root}/*.{ext}')
        if not cbfs:
            continue
        allFilesR = run(root,ponifile, stdevs, maskfile, scale, threshold, polarisation, gainFile, nbins, ext,outdir, save2d=save2d, 
            saveMasks=saveMasks, cpp= cpp, allFiles=allFilesR, averaging=averaging)
        
    return allFilesR

if __name__ == '__main__':
    pass
    #run(datadir, poni,  stdevs, maskfile, scale, threshold, polarisation,gainFile, save2d=save2d, saveMasks=saveMasks, nbins=nbins)