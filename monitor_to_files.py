import fabio
import os
from glob import glob
import pandas as pd
'''
for putting monitor count values from the log files into the image file headers
'''
direc = r'C:\Users\kenneth1a\Documents\beamlineData\a311189_pdf\01_RbCuF3_Q0p5\pdf_800C\xrd'
os.chdir(direc)
def addExtension(fname):
    if not fname.endswith('.cbf'):
        fname +=  '.cbf'
    return fname

monitorFile = glob('*.dat')[0] #'empty_B1mm_pdf_001_dtx_0_counts.dat'
df = pd.read_csv(monitorFile,sep = ' ')
df['file'] = df['file'].apply(addExtension)

cbfs = glob('*.cbf')
for file in cbfs:
    data = fabio.open(file)
    monitorCounts = df[df['file'] == file]['monitor'].values[0]
    if 'Flux' in data.header["_array_data.header_contents"]:
        continue
    data.header["_array_data.header_contents"] += f'\r\n# Flux {monitorCounts}'
    data.save(file)
    