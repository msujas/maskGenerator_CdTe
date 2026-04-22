import maskGenerator_gainMaps
from glob import glob
posList = [f'pos{i:02d}' for i in range(1,12)]
gainFile = None #r'Z:\visitor\ch6987\bm31\20240326\pylatus/gainMap_filtered_kpm_2024-03-28.edf' #None
ponidir = fr'C:\Users\kenneth1a\Documents\beamlineData\Dec2025_gain_C60\gainmap90'
for c,pos in enumerate(posList,1):

    direc = fr'{ponidir}/glassRod/{pos}'
    dest = direc
    mask = fr'{ponidir}/thlimMask/{pos}_mask.edf'
    poni = f'{ponidir}/{pos}_MD.poni'
    print(direc)
    maskGenerator_gainMaps.run(direc,mask,gainFile, poni = poni,avdir = 'average',dest = dest, doMonitor=True)
    