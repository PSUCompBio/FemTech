import sys
import json
import linecache
import numpy as np

import csv

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

if (not 'compute-injury-criteria' in inputJson['simulation']) or inputJson['simulation']['compute-injury-criteria']:
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

if 'output-nodes' in inputJson['simulation'] or 'output-elements' in inputJson['simulation']:
    # Open plot column file 
    plotFile = 'plot_' + uid + '.dat'
    with open(plotFile, "r") as csv_file:
        csv_reader = csv.DictReader(csv_file, delimiter=' ')
        # Get count to pre-allocate arrays
        nRow = sum(1 for _ in csv_reader)
        # Read the time arrays
        time = np.zeros(nRow)
        csv_file.seek(0)
        # Go to next item to skip header
        next(csv_reader)
        for n, row in enumerate(csv_reader):
            time[n] = row['#Time']
        outputJson['plot'] = {}
        outputJson['plot']['time'] = time.tolist()
        # Read all the nodes
        xVal = np.zeros(nRow)
        yVal = np.zeros(nRow)
        zVal = np.zeros(nRow)
        for n, node in enumerate(inputJson['simulation']['output-nodes']):
            nodeStr = 'Node'+'%08d'%node+'-Disp'
            csv_file.seek(0)
            # Go to next item to skip header
            next(csv_reader)
            for i, row in enumerate(csv_reader):
                xVal[i] = row[nodeStr+'X']
                yVal[i] = row[nodeStr+'Y']
                zVal[i] = row[nodeStr+'Z']
            outputJson['plot']['nodal-displacement-'+str(n)] = {}
            outputJson['plot']['nodal-displacement-'+str(n)]['x'] = xVal.tolist();
            outputJson['plot']['nodal-displacement-'+str(n)]['y'] = yVal.tolist();
            outputJson['plot']['nodal-displacement-'+str(n)]['z'] = zVal.tolist();
        for n, element in enumerate(inputJson['simulation']['output-elements']):
            elemStr = 'Elem'+'%08d'%element+'-Stre'
            csv_file.seek(0)
            # Go to next item to skip header
            next(csv_reader)
            for i, row in enumerate(csv_reader):
                xVal[i] = row[elemStr+'P']
                yVal[i] = row[elemStr+'S']
            outputJson['plot']['element-stress-'+str(n)] = {}
            outputJson['plot']['element-stress-'+str(n)]['principal'] = xVal.tolist();
            outputJson['plot']['element-stress-'+str(n)]['shear'] = yVal.tolist();

# print(outputJson)
# jstr = json.dumps(outputJson, indent = 2, sort_keys=False)
# print(jstr)
# jsonOutputFile = uid + '_output_modified.json'
with open(jsonOutputFile, 'w') as outfile:
    json.dump(outputJson, outfile, indent = 2, sort_keys=False)
