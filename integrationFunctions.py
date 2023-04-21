import fabio
import numpy as np
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

