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
   from integrationFunctions import clearPyFAI_header, gainCorrection 
else:
    from . import clearPyFAI_header, gainCorrection

poni = r'W:\simonWintersteller/Si_00_0tilt.poni'
polarisation = 0.99
stdevs = 4
datadir = r'W:\simonWintersteller\pdf_0_0tilt_33thresh/'
maskfile = r'W:\simonWintersteller/mask_0_00tilt2024_4rows.edf'
scale = 10**9
gainFile = r'W:\simonWintersteller/gainMap_filtered_kpm_2025-02-05.edf'


def generateMask(dataarray, binarray,nbins, stdevs, basemask):
    array = np.empty(shape = dataarray.shape)
    newmask = basemask.copy()
    for n in range(nbins):
        wherebin = np.where(binarray == n)
        array =  dataarray[wherebin]
        maskBin = basemask[wherebin]
        stdev = np.nanstd(array)
        median = np.nanmedian(array)
        maskBin = np.where(array > median + stdevs*stdev, 1, maskBin)
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

def run(datadir, ponifile, polarisation, stdevs, maskfile, scale, gainFile = None, nbins = 800,outdir = 'xye'):
    os.chdir(datadir)
    if not os.path.exists(f'{datadir}/{outdir}/'):
        os.makedirs(f'{datadir}/{outdir}/')
    basemask = fabio.open(maskfile).data
    cbfs = glob(f'*.cbf')
    array2th, polarray, saArray, binarray = generateDetArrays(cbfs[0],ponifile,nbins)
    poni = pyFAI.load(ponifile)
    maskdir = f'{datadir}/masks{stdevs}/'
    if gainFile != None:
        gainMap = fabio.open(gainFile).data
    os.makedirs(maskdir,exist_ok=True)
    for file in cbfs:
        print(file)
        dataarray = makeArray(file,polarray, saArray,gainMap)
        xyefile = file.replace('.cbf','.xye')
        outfile = f'{outdir}/{xyefile}'
        header = CbfHeader(file)
        monitorCounts = header['Flux']
        mask = generateMask(dataarray,binarray,nbins, stdevs, basemask)
        maskim = fabio.edfimage.EdfImage(mask)
        maskim.save(f'{maskdir}/{file}'.replace('.cbf','.edf'))
        normArray = (dataarray/monitorCounts) * scale
        poni.integrate1d(normArray,mask = mask, filename = outfile, polarization_factor = polarisation,
                        unit = '2th_deg', correctSolidAngle = True, method = 'bbox', npt = 5000, 
                        error_model = 'poisson', safe = False)
        clearPyFAI_header(outfile)

if __name__ == '__main__':
    run(datadir, poni, polarisation, stdevs, maskfile, scale,gainFile)