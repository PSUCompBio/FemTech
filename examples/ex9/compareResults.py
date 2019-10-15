import os
import sys
import numpy as np

def getLastLine(fname):
    with open(fname, 'rb') as f:
        f.seek(-2, os.SEEK_END)
        while f.read(1) != b'\n':
            f.seek(-2, os.SEEK_CUR)
        lastLine = f.readline().decode()
        array = map(float, lastLine.split())
        return np.array(array)

femTech = getLastLine('plot.dat')
abaqus = getLastLine('abaqus/abaqus.rpt')

tF = femTech[0]; xF = femTech[1:];
tA = abaqus[0]; xA = abaqus[10:13];

dT = np.abs((tF-tA)*100.0/tA);
dX = np.abs((xF-xA)*100.0/xA);

if (dT < 1.0):
    dMax = dX.max()
    if (dMax > 3.0):
        raise sys.exit(1)
else:
    raise sys.exit(1)
