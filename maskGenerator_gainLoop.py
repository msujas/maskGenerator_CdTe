import maskGenerator_gainMaps
from glob import glob
posList = [f'pos{i:02d}' for i in range(1,12)]
gainFile = r'Z:\visitor\ch6987\bm31\20240326\pylatus/gainMap_filtered_kpm_2024-03-28.edf' #None
ponidir = fr'Z:\visitor\ch6987\bm31\20240326\pylatus\gainmap'
for c,pos in enumerate(posList,1):

    direc = fr'Z:\visitor\ch6987\bm31\20240326\pylatus\gainmap/glassRod/{pos}'
    dest = direc.replace(r'X:\users\a311217',r'Z:\visitor\a311217\bm31\20240129\pylatus')
    mask = fr'{ponidir}/{pos}_mask.edf'
    poni = f'{ponidir}/{pos}_MD.poni'
    print(direc)
    maskGenerator_gainMaps.run(direc,mask,gainFile, poni = poni,avdir = 'average',dest = dest, doMonitor=True)
    