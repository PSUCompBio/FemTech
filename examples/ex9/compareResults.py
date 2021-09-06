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
abaqus = getLastLine('ls-dyna/ogdenDisplacement.txt')

print(femTech)
print(abaqus)

tF = femTech[0]; xF = femTech[1:];
tA = abaqus[0]; xA = abaqus[1:];

dT = np.abs((tF-tA)*100.0/tA);
dX = np.abs((xF-xA)*100.0/xA);

if (dT < 1.0):
    dMax = dX.max()
    if (dMax > 5.0):
        raise sys.exit(1)
else:
    raise sys.exit(1)
print "Validation Sucessful : Ex9"
