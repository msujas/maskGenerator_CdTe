import maskGenerator_gainMaps
from glob import glob
posList = [f'pos{i:02d}' for i in range(1,12)]
gainFile = None#r'Z:\visitor\a311217\bm31\20240129\pylatus\gainmap/gainMap_thr30keV_filtered_kpm_2024-02-01.edf'
ponidir = fr'Z:\visitor\a311217\bm31\20240129\pylatus\gainmap'
for c,pos in enumerate(posList,1):

    direc = fr'X:\users\a311217\gainmap\glassRod/{pos}'
    dest = direc.replace(r'X:\users\a311217',r'Z:\visitor\a311217\bm31\20240129\pylatus')
    mask = fr'{ponidir}/thlimMask/{pos}_mask.edf'
    poni = f'{ponidir}/{pos}_MD.poni'
    print(direc)
    maskGenerator_gainMaps.run(direc,mask,gainFile, poni = poni,avdir = 'average',dest = dest, doMonitor=True)
    