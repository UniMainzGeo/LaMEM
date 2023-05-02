import sys
import os

flist = os.listdir('.')

for file in flist:

    ext = file.split(".")[-1]

    if ext == 'cpp':

        cmd = 'grep -v "__FUNCT__" '+ file + ' > temp && mv temp ' + file

        # print(cmd)

        os.system(cmd)
