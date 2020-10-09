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

if 'xt' in linearAcc:
    maxG = np.max(np.sqrt(np.array(linearAcc['xv'])**2+np.array(linearAcc['yv'])**2+np.array(linearAcc['zv'])**2))/9.81
    maxT = np.max(np.sqrt(np.array(angularAcc['xv'])**2+np.array(angularAcc['yv'])**2+np.array(angularAcc['zv'])**2))
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
