from cryio.cbfimage import CbfHeader
import os
from glob import glob
#import shutil

sourcedir = r'X:\users\a311217'
destdir = r'Z:\visitor\a311217\bm31\20240129\pylatus'

def run(sourcedir,destdir):
    for root,dirs, files in os.walk(sourcedir):
        cbffiles = glob(f'{root}/*.cbf')
        if len(cbffiles) == 0:
            continue
        
        destfiles = [file.replace(sourcedir,destdir) for file in cbffiles]
        for file,destfile in zip(cbffiles,destfiles):
            header = CbfHeader(file)
            if not 'Flux' in header:
                continue
            flux = int(header['Flux'])
            print(file)
            print(destfile)
            if not os.path.exists(destfile):
                #shutil.copy(file, destfile)
                continue
            destheader = CbfHeader(destfile)
            destheader.header['Flux'] = flux
            destheader.save_cbf(destfile)
            
if __name__ == '__main__':
    run(sourcedir,destdir) 
    