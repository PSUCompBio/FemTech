import sys
import json
import linecache
import numpy as np

jsonInputFile = sys.argv[1]
inputJson = json.loads(open(jsonInputFile).read())
uid = inputJson["uid"]
jsonOutputFile = uid + '_output.json'
linearAcc = inputJson['simulation']['linear-acceleration']
angularAcc = inputJson['simulation']['angular-acceleration']

# Check if time array is present, common for all accelerations 
if 'time-all' in inputJson['simulation']:
    maxG = np.max(np.sqrt(np.array(linearAcc['xv'])**2+np.array(linearAcc['yv'])**2+np.array(linearAcc['zv'])**2))/9.81
    maxT = np.max(np.sqrt(np.array(angularAcc['xv'])**2+np.array(angularAcc['yv'])**2+np.array(angularAcc['zv'])**2))
elif 'xt' in linearAcc: # If all sensor values have corresponding time arrays
    # If they are equal compute maximum directly, else interpolate and compute
    if (np.array_equal(np.array(linearAcc['xt']), np.array(linearAcc['yt'])) 
            and np.array_equal(np.array(linearAcc['xt']), np.array(linearAcc['zt']))):
        maxG = np.max(np.sqrt(np.array(linearAcc['xv'])**2+np.array(linearAcc['yv'])**2+np.array(linearAcc['zv'])**2))/9.81
    else:
        nCount = 5*max([len(linearAcc['xt']), len(linearAcc['yt']), len(linearAcc['zt'])])
        tmin = min(linearAcc['xt'][0], linearAcc['yt'][0], linearAcc['zt'][0])
        tmax = min(linearAcc['xt'][-1], linearAcc['yt'][-1], linearAcc['zt'][-1])
        tInter = np.linspace(tmin, tmax, nCount)
        xV = np.interp(tInter, linearAcc['xt'], linearAcc['xv'])
        yV = np.interp(tInter, linearAcc['yt'], linearAcc['yv'])
        zV = np.interp(tInter, linearAcc['zt'], linearAcc['zv'])
        maxG = np.max(np.sqrt(xV**2+yV**2+zV**2))/9.81
    if (np.array_equal(np.array(angularAcc['xt']), np.array(angularAcc['yt']))
            and np.array_equal(np.array(angularAcc['xt']), np.array(angularAcc['zt']))):
        maxT = np.max(np.sqrt(np.array(angularAcc['xv'])**2+np.array(angularAcc['yv'])**2+np.array(angularAcc['zv'])**2))
    else:
        nCount = 5*max([len(angularAcc['xt']), len(angularAcc['yt']), len(angularAcc['zt'])])
        tmin = min(angularAcc['xt'][0], angularAcc['yt'][0], angularAcc['zt'][0])
        tmax = min(angularAcc['xt'][-1], angularAcc['yt'][-1], angularAcc['zt'][-1])
        tInter = np.linspace(tmin, tmax, nCount)
        xV = np.interp(tInter, angularAcc['xt'], angularAcc['xv'])
        yV = np.interp(tInter, angularAcc['yt'], angularAcc['yv'])
        zV = np.interp(tInter, angularAcc['zt'], angularAcc['zv'])
        maxT = np.max(np.sqrt(xV**2+yV**2+zV**2))
else:
    maxG = np.sqrt(np.array(linearAcc[0])**2+np.array(linearAcc[1])**2+np.array(linearAcc[2])**2)
    if isinstance(angularAcc, list):
        maxT = np.sqrt(np.array(angularAcc[0])**2+np.array(angularAcc[1])**2+np.array(angularAcc[2])**2)
    else:
        maxT = np.fabs(angularAcc)

outputJson = json.loads(open(jsonOutputFile).read())
meshType = outputJson["output-file"].split('_')[0]
cellDataFile = meshType + '_cellcentres.txt'

# Part list 
partList = ["msc", "stem", "cerebellum", "frontal", "parietal", "occipital",
        "temporal"]

# Populate region from cell centres file
maxInjuryMetrics = ["principal-max-strain", "principal-min-strain",
        "maximum-shear-strain", "maximum-PSxSR"];
threshInjuryMetrics = ["CSDM-10", "CSDM-15", "CSDM-30", "CSDM-5", 
        "MPS-95"];

# Loop over metric with single element output
for metric in maxInjuryMetrics:
    maxLocation = int(outputJson[metric]["global-element-id"])
    maxLine = linecache.getline(cellDataFile, maxLocation).split()
    outputJson[metric]["brain-region"] = maxLine[4]
    # overwrite location with femtech-reference co-ordinate
    outputJson[metric]["location"] = \
            [float(maxLine[1]), float(maxLine[2]), float(maxLine[3])]
    outputJson[metric].pop("global-element-id", None)

# Loop over metric with multiple element output
for metric in threshInjuryMetrics:
    maxLocation = outputJson[metric]["global-element-id"]
    for part in partList:
        outputJson[metric][part] = []
    for loc in maxLocation:
        maxLine = linecache.getline(cellDataFile, int(loc)).split()
        outputJson[metric][maxLine[4]].append([float(maxLine[1]), \
                float(maxLine[2]), float(maxLine[3])])
    outputJson[metric].pop("global-element-id", None)

# Add max quantities 
outputJson["max-linear-acc-g"] = maxG
outputJson["max-angular-acc-rads2"] = maxT

# jstr = json.dumps(outputJson, indent = 2, sort_keys=False)
# print(jstr)
with open(jsonOutputFile, 'w') as outfile:
    json.dump(outputJson, outfile, indent = 2, sort_keys=False)
