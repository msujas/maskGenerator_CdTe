import os, shutil

sourcedir = r'X:\users\a311217'
destdir = r'Z:\visitor\a311217\bm31\20240129\pylatus'
for root,dirs,files in os.walk(sourcedir):
    sourcefiles = [f'{root}/{file}' for file in files]
    for file in sourcefiles:
        destfile = file.replace(sourcedir,destdir)
        destsubdir = os.path.dirname(destfile)
        if not os.path.exists(destfile):
            if not os.path.exists(destsubdir):
                os.makedirs(destsubdir)
            shutil.copyfile(file,destfile)
            print(destfile)

