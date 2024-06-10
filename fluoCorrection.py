import numpy as np
import pyFAI
import matplotlib.pyplot as plt

def solidAngle(poni1,poni2, d, px, py):
    psize = 172e-6 #m
    xpos = px*psize
    ypos = py*psize
    angle1 = np.arctan((np.abs(ypos-poni1)+psize/2)/d) - np.arctan((np.abs(ypos-poni1)-psize/2)/d)
    angle2 = np.arctan((np.abs(xpos-poni2)+psize/2)/d) - np.arctan((np.abs(xpos-poni2)-psize/2)/d)
    return angle1*angle2

det = np.empty(shape = (1679,1475,2))

for y in range(len(det)):
    for x in range(len(det[0])):
        det[y,x] = [y,x]
det = det.astype('uint16')

def solidAngleMap(poni1,poni2,d):
    return solidAngle(poni1,poni2,d,det[:,:,1],det[:,:,0])

poni1 = 0.23629
poni2 = 0.12924
d = 0.186

#saMap = solidAngleMap(poni1,poni2, d)

def fluoCorrection(poni1,poni2,d, fluoK):
    return fluoK*solidAngleMap(poni1,poni2,d)

fCorr=  fluoCorrection(poni1,poni2,d,10**7)

plt.imshow(fCorr)
plt.colorbar()
plt.show()



