import pandas as pd
import json, re
import matplotlib.pyplot as plt
import numpy as np
import sys

# Set plot parameters
params = {'legend.fontsize': 'xx-large',
         'axes.labelsize': 'x-large',
         'axes.titlesize':'xx-large',
         'xtick.labelsize':'x-large',
         'ytick.labelsize':'x-large'}
plt.rcParams.update(params)

fileName = sys.argv[1]
# Remove comments in json string before parsing
config = json.loads(re.sub("//.*","", open(fileName).read(),flags=re.MULTILINE))

uid = config["uid"]
# Read the main output file
trace = pd.read_csv('./plot_'+uid+'.dat', delimiter=' ')
logFile = './femtech_'+uid+'.log'

settings = config["simulation"]
outputNodes = np.array(settings["output-nodes"])
C = np.char.mod('%d', outputNodes)
nC = outputNodes.size
coord = np.zeros((nC, 3), dtype=float)

# Read co-ordinates from log file
with open(logFile, 'r') as f:
    startPos = f.tell()
    for index, cElem in enumerate(C):
        lineSearch = ': Node ID : '+cElem 
        for line in f:
            if lineSearch in line:
                cordString = re.search('\((.*)\)', line)
                coord[index, :] = np.array(cordString.group(1).split(', '), np.double)
                break
        f.seek(startPos)

C = np.char.mod('%08d', outputNodes)

xLimMax = 75.0
yLimMax = 50.0
zLimMax = 100.0
xticks = np.arange(-xLimMax, xLimMax+1, 25.0)
yticks = np.arange(-yLimMax, yLimMax+1, 25.0)
zticks = np.arange(zLimMax, -1, -25.0)

################################################################
# Sagittal plot
plt.figure(figsize = (12, 12*zLimMax/(2*xLimMax)))
plt.grid(color='k', linestyle='-', linewidth=0.05)
axis = plt.gca()
axis.set(xticks = xticks, xticklabels = xticks)
axis.set(yticks = zticks, yticklabels = zticks)
axis.xaxis.set_ticks_position('top')
axis.xaxis.set_label_position('top')
# Move left y-axis and bottim x-axis to centre, passing through (0,0)
axis.spines['left'].set_position('center')
# Eliminate upper and right axes
axis.spines['right'].set_color('none')
axis.spines['bottom'].set_color('none')
# Show ticks in the left and lower axes only
axis.xaxis.set_ticks_position('top')
axis.yaxis.set_ticks_position('left')
axis.yaxis.set_label_coords(-0.02, 0.5)
axis.axis('equal')
axis.set(xlim=(-xLimMax-5, xLimMax+5), ylim=(zLimMax+10, 0.0))

for index, cElem in enumerate(C):
    cName = 'Node'+cElem+'-Disp'
    X = (np.array(trace[cName+'X'])+coord[index, 0])*1000
    Z = (np.array(trace[cName+'Z'])+coord[index, 2])*1000
    plt.plot(X, Z, zorder=10)
    plt.text(coord[index, 0]*1000, (coord[index, 2])*1000-2.0, 'C'+str(index+1),
            fontsize=8, zorder=1, color='silver')

# plt.xlim(-xLimMax-5, xLimMax+5)
# plt.ylim(zLimMax, 0.0)
plt.xlabel('X (mm)')
plt.ylabel('Z (mm)')
plt.tight_layout()
plt.savefig('C288T3_sagittal.png')
plt.close()

################################################################
# Coronal plot
plt.figure(figsize = (12, 12*zLimMax/(2*yLimMax)))
plt.grid(color='k', linestyle='-', linewidth=0.05)
axis = plt.gca()
axis.set(xticks = yticks, xticklabels = yticks)
axis.set(yticks = zticks, yticklabels = zticks)
axis.xaxis.set_ticks_position('top')
axis.xaxis.set_label_position('top')
# Move left y-axis and bottim x-axis to centre, passing through (0,0)
axis.spines['left'].set_position('center')
# Eliminate upper and right axes
axis.spines['right'].set_color('none')
axis.spines['bottom'].set_color('none')
# Show ticks in the left and lower axes only
axis.xaxis.set_ticks_position('top')
axis.yaxis.set_ticks_position('left')
axis.yaxis.set_label_coords(-0.02, 0.5)
axis.axis('equal')
axis.set(xlim=(-yLimMax-5, yLimMax+5), ylim=(zLimMax+10, 0.0))

for index, cElem in enumerate(C):
    cName = 'Node'+cElem+'-Disp'
    Y = (np.array(trace[cName+'Y'])+coord[index, 1])*1000
    Z = (np.array(trace[cName+'Z'])+coord[index, 2])*1000
    plt.plot(Y, Z, zorder=10)
    plt.text(coord[index, 1]*1000, (coord[index, 2])*1000-2.0, 'C'+str(index+1),
            fontsize=8, zorder=1, color='silver')

plt.xlabel('Y (mm)')
plt.ylabel('Z (mm)')
plt.tight_layout()
plt.savefig('C288T3_coronal.png')
plt.close()

################################################################
# Horizontal plot
plt.figure(figsize = (12, 12*yLimMax/xLimMax))
plt.grid(color='k', linestyle='-', linewidth=0.05)
axis = plt.gca()
axis.set(xticks = xticks, xticklabels = xticks)
axis.set(yticks = yticks, yticklabels = yticks)
axis.xaxis.set_ticks_position('top')
axis.xaxis.set_label_position('top')
# Move left y-axis and bottim x-axis to centre, passing through (0,0)
axis.spines['left'].set_position('center')
axis.spines['top'].set_position('center')
# Eliminate upper and right axes
axis.spines['right'].set_color('none')
axis.spines['bottom'].set_color('none')
# Show ticks in the left and lower axes only
axis.xaxis.set_ticks_position('top')
axis.yaxis.set_ticks_position('left')
axis.yaxis.set_label_coords(-0.02, 0.5)
axis.xaxis.set_label_coords(0.5, 1.02)
axis.axis('equal')
axis.set(xlim=(-xLimMax-5, xLimMax+5), ylim=(-yLimMax-5, yLimMax+5))

for index, cElem in enumerate(C):
    cName = 'Node'+cElem+'-Disp'
    X = (np.array(trace[cName+'X'])+coord[index, 0])*1000
    Y = (np.array(trace[cName+'Y'])+coord[index, 1])*1000
    plt.plot(X, Y, zorder=10)
    plt.text(coord[index, 0]*1000, (coord[index, 1])*1000-2.0, 'C'+str(index+1),
            fontsize=8, zorder=1, color='silver')

plt.xlabel('X (mm)')
plt.ylabel('Y (mm)')
plt.tight_layout()
plt.savefig('C288T3_horizontal.png')
plt.close()

################################################################
# Relative displacment plots
################################################################
T = np.array(trace['#Time'])*1000
fig, axs = plt.subplots(3, 7, figsize=(21,9))
for index, cElem in enumerate(C[0:7]):
    cName = 'Node'+cElem+'-Disp'
    X = np.array(trace[cName+'X'])*1000
    Z = np.array(trace[cName+'Y'])*1000
    Y = np.array(trace[cName+'Y'])*1000
    axs[0, index].plot(T, X)
    title = 'NDTx'+str(index+1)
    axs[0, index].set_title(title)
    axs[0, index].set_ylabel('X-Disp(mm)')
    axs[1, index].plot(T, Y)
    axs[1, index].set_title(title)
    axs[1, index].set_ylabel('Y-Disp(mm)')
    axs[2, index].plot(T, Z)
    axs[2, index].set_title(title)
    axs[2, index].set_ylabel('Z-Disp(mm)')

# Set global plot properties
xticks = np.arange(0, 40+1, 20)
yticks = np.arange(-5, 5+1, 5)
for axsE in axs.flat:
    axsE.set(xlim=(0, 40), ylim=(-5, 5))
    axsE.set(xticks = xticks, xticklabels = xticks)
    axsE.set(yticks = yticks, yticklabels = yticks)
    axsE.set_xlabel('Time(ms)')

plt.tight_layout()
plt.savefig('C288T3_relativeDisplacement_1.png')
plt.close()

if nC > 7:
    fig, axs = plt.subplots(3, 7, figsize=(21,9))
    for index, cElem in enumerate(C[7:14]):
        cName = 'Node'+cElem+'-Disp'
        X = np.array(trace[cName+'X'])*1000
        Z = np.array(trace[cName+'Y'])*1000
        Y = np.array(trace[cName+'Y'])*1000
        axs[0, index].plot(T, X)
        title = 'NDTx'+str(index+8)
        axs[0, index].set_title(title)
        axs[0, index].set_ylabel('X-Disp(mm)')
        axs[1, index].plot(T, Y)
        axs[1, index].set_title(title)
        axs[1, index].set_ylabel('Y-Disp(mm)')
        axs[2, index].plot(T, Z)
        axs[2, index].set_title(title)
        axs[2, index].set_ylabel('Z-Disp(mm)')

    # Set global plot properties
    xticks = np.arange(0, 40+1, 20)
    yticks = np.arange(-5, 5+1, 5)
    for axsE in axs.flat:
        axsE.set(xlim=(0, 40), ylim=(-5, 5))
        axsE.set(xticks = xticks, xticklabels = xticks)
        axsE.set(yticks = yticks, yticklabels = yticks)
        axsE.set_xlabel('Time(ms)')

    plt.tight_layout()
    plt.savefig('C288T3_relativeDisplacement_2.png')
    plt.close()

################################################################
# Plot principal strains
################################################################
outputElems = np.array(settings["output-elements"])
C = np.char.mod('%08d', outputElems)
nC = outputElems.size

xLimMax = 120
yLimMax = 0.2
xticks1 = np.arange(0, xLimMax+1, 30)
xticks2 = np.arange(0, 40+1, 10)
yticks = np.arange(0, yLimMax+0.01, 0.05)

for index, cElem in enumerate(C):
    plt.figure()
    cName = 'Elem'+cElem+'-Stre'
    X1 = np.array(trace[cName+'P'])
    X2 = np.array(trace[cName+'S'])
    plt.plot(T, X1, linewidth=2.0, label='Principal')
    plt.plot(T, X2, linewidth=2.0, label='Shear')
    plt.legend(loc='best')
    plt.xlabel('Time (ms)')
    plt.ylabel('Principal Strain')
    plt.gca().set(xlim=(0, 120), ylim=(0, 0.2))
    plt.gca().set(xticks = xticks1, xticklabels = xticks1)
    plt.gca().set(yticks = yticks)
    plt.tight_layout()
    plt.savefig('C288T3_C'+str(index+1)+'_Strain.png')
    plt.close()

    plt.figure()
    plt.plot(T, X1, linewidth=2.0, label='Principal')
    plt.plot(T, X2, linewidth=2.0, label='Shear')
    plt.legend(loc='best')
    plt.xlabel('Time (ms)')
    plt.ylabel('Principal Strain')
    plt.gca().set(xlim=(0, 40), ylim=(0, 0.2))
    plt.gca().set(xticks = xticks2, xticklabels = xticks2)
    plt.gca().set(yticks = yticks)
    plt.tight_layout()
    plt.savefig('C288T3_C'+str(index+1)+'_Strain_Small.png')
    plt.close()
