import sys
import json
import linecache
import numpy as np

import csv

# Part list 
partList = ["msc", "stem", "cerebellum", "frontal", "parietal", "occipital",
        "temporal"]

# Code assumes maximum-time is set to the last value of the time array
jsonInputFile = sys.argv[1]
inputJson = json.loads(open(jsonInputFile).read())
uid = inputJson["event_id"]
jsonOutputFile = uid + '_output.json'

computeAngAcc = False
if 'angular-acceleration' in inputJson['simulation']:
    computeAngAcc = True
    angularAcc = inputJson['simulation']['angular-acceleration']

computeAngVel = False
if 'angular-velocity' in inputJson['simulation']:
    computeAngVel = True
    angularVel = inputJson['simulation']['angular-velocity']

pressureSimulation = False
if 'pressure' in inputJson['simulation']:
    pressureSimulation = True

if not pressureSimulation:
    linearAcc = inputJson['simulation']['linear-acceleration']
    # Check if time array is present, common for all accelerations 
    if 'time-all' in inputJson['simulation']:
        maxG = np.max(np.sqrt(np.array(linearAcc['xv'])**2+np.array(linearAcc['yv'])**2+np.array(linearAcc['zv'])**2))/9.81
        if computeAngAcc:
            maxT = np.max(np.sqrt(np.array(angularAcc['xv'])**2+np.array(angularAcc['yv'])**2+np.array(angularAcc['zv'])**2))
        if computeAngVel:
            maxAV = np.max(np.sqrt(np.array(angularVel['xv'])**2+np.array(angularVel['yv'])**2+np.array(angularVel['zv'])**2))
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
        if computeAngAcc:
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

        if computeAngVel:
            if (np.array_equal(np.array(angularVel['xt']), np.array(angularVel['yt']))
                    and np.array_equal(np.array(angularVel['xt']), np.array(angularVel['zt']))):
                maxAV = np.max(np.sqrt(np.array(angularVel['xv'])**2+np.array(angularVel['yv'])**2+np.array(angularVel['zv'])**2))
            else:
                nCount = 5*max([len(angularVel['xt']), len(angularVel['yt']), len(angularVel['zt'])])
                tmin = min(angularVel['xt'][0], angularVel['yt'][0], angularVel['zt'][0])
                tmax = min(angularVel['xt'][-1], angularVel['yt'][-1], angularVel['zt'][-1])
                tInter = np.linspace(tmin, tmax, nCount)
                xV = np.interp(tInter, angularVel['xt'], angularVel['xv'])
                yV = np.interp(tInter, angularVel['yt'], angularVel['yv'])
                zV = np.interp(tInter, angularVel['zt'], angularVel['zv'])
                maxAV = np.max(np.sqrt(xV**2+yV**2+zV**2))
    else:
        maxG = np.sqrt(np.array(linearAcc[0])**2+np.array(linearAcc[1])**2+np.array(linearAcc[2])**2)
        if computeAngAcc:
            if isinstance(angularAcc, list):
                maxT = np.sqrt(np.array(angularAcc[0])**2+np.array(angularAcc[1])**2+np.array(angularAcc[2])**2)
            else:
                maxT = np.fabs(angularAcc)
else: # pressure simulation
    pressureJson = inputJson['simulation']['pressure']
    maxP = np.max(np.abs(np.array(pressureJson['head'])+np.array(pressureJson['shoulder'])+np.array(pressureJson['chest'])))/3.0

outputJson = json.loads(open(jsonOutputFile).read())

if (not 'compute-injury-criteria' in inputJson['simulation']) or inputJson['simulation']['compute-injury-criteria']:
    meshFileParts = outputJson["output-file"].split('_')
    if 'coarse' in meshFileParts:
        meshType = 'coarse'
    else:
        meshType = 'fine'
    if 'female' in meshFileParts:
        meshSuffix = 'female'
    else:
        meshSuffix = 'male'
    cellDataFile = meshType + '_cellcentresandvol_'+ meshSuffix +'.txt'

    # Populate region from cell centres file
    maxInjuryMetrics = ["principal-max-strain", "principal-min-strain",
            "maximum-shear-strain", "maximum-PSxSR", "principal-max-pressure"];
    threshInjuryMetrics = ["CSDM-10", "CSDM-15", "CSDM-30", "CSDM-5", "CSDM-3"];

    # Loop over metric with single element output
    for metric in maxInjuryMetrics:
        maxLocation = int(outputJson[metric]["global-element-id"])
        maxLine = linecache.getline(cellDataFile, maxLocation).split()
        outputJson[metric]["brain-region"] = maxLine[4]
        # overwrite location with femtech-reference co-ordinate
        outputJson[metric]["location"] = \
                [float('%.3g' % float(maxLine[1])), float('%.3g' % float(maxLine[2])), \
                            float('%.3g' % float(maxLine[3]))]
        outputJson[metric].pop("global-element-id", None)

    # Loop over metric with multiple element output
    for metric in threshInjuryMetrics:
        jsonOutputFileMetric = metric + '.json'
        outputJsonMetric = json.loads(open(jsonOutputFileMetric).read())
        maxLocation = outputJsonMetric["global-element-id"]
        for part in partList:
            outputJsonMetric[part] = []
        for loc in maxLocation:
            maxLine = linecache.getline(cellDataFile, int(loc)).split()
            outputJsonMetric[maxLine[4]].append([float('%.3g' % float(maxLine[1])), \
                    float('%.3g' % float(maxLine[2])), float('%.3g' % float(maxLine[3]))])
        outputJsonMetric.pop("global-element-id", None)
        with open(jsonOutputFileMetric, 'w') as outfile:
            json.dump(outputJsonMetric, outfile, indent = 2, sort_keys=False)

# Add max quantities 
if not pressureSimulation:
    outputJson["max-linear-acc-g"] = float('%.3g' % maxG)
    if computeAngAcc:
        outputJson["max-angular-acc-rads2"] = float('%.3g' % maxT)
    if computeAngVel:
        outputJson["max-angular-vel-rads"] = float('%.3g' % maxAV)
else:
    outputJson["max-pressure"] = float('%.3g' % maxP)

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
        if 'output-nodes' in inputJson['simulation']:
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
        if 'output-elements' in inputJson['simulation']:
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

if 'impact-time' in inputJson:
    outputJson['impact-time'] = inputJson['impact-time']
if 'impact-date' in inputJson:
    outputJson['impact-date'] = inputJson['impact-date']

# Include region wise MPS values to output.json
if (not 'compute-injury-criteria' in inputJson['simulation']) or inputJson['simulation']['compute-injury-criteria']:
    outputJson["region-wise-mps"] = {}
    # Set all region MPS value to zero
    for part in partList:
        outputJson["region-wise-mps"][part] = 0.0

    # Create and write MPS-95.json
    outputJsonMetric = {}
    outputJsonMetric["time"] = outputJson["MPS-95"]["time"]
    outputJsonMetric["value"] = outputJson["MPS-95"]["value"]
    for part in partList:
        outputJsonMetric[part] = []

    with open('MPSfile.dat') as mpsFile, open(cellDataFile) as regionFile:
        for mpsLine, cellDataLine in zip(mpsFile, regionFile):
            mpsLineSplit = cellDataLine.split();
            part = mpsLineSplit[4]
            if part != 'skull' and part != 'csf':
                mpsValue = float(mpsLine.split(',')[1])
                if mpsValue > outputJson["region-wise-mps"][part]:
                    outputJson["region-wise-mps"][part] = float('%.3g' % mpsValue)
                if mpsValue >= outputJsonMetric["value"]:
                    outputJsonMetric[part].append([float('%.3g' % float(mpsLineSplit[1])), \
                    float('%.3g' % float(mpsLineSplit[2])), float('%.3g' % float(mpsLineSplit[3]))])
    jsonOutputFileMPS95 = 'MPS-95.json'
    with open(jsonOutputFileMPS95, 'w') as outfile:
        json.dump(outputJsonMetric, outfile, indent = 2, sort_keys=False)

# Include region wise MPr values to output.json
if (not 'compute-injury-criteria' in inputJson['simulation']) or inputJson['simulation']['compute-injury-criteria']:
    outputJson["region-wise-mpr"] = {}
    # Set all region MPr value to zero
    for part in partList:
        outputJson["region-wise-mpr"][part] = 0.0

    # Create and write MPr-95.json
    outputJsonMetric = {}
    outputJsonMetric["time"] = outputJson["MPr-95"]["time"]
    outputJsonMetric["value"] = outputJson["MPr-95"]["value"]
    for part in partList:
        outputJsonMetric[part] = []

    with open('MPrFile.dat') as mprFile, open(cellDataFile) as regionFile:
        for mprLine, cellDataLine in zip(mprFile, regionFile):
            mprLineSplit = cellDataLine.split();
            part = mprLineSplit[4]
            if part != 'skull' and part != 'csf':
                mprValue = float(mprLine.split(',')[1])
                if mprValue > outputJson["region-wise-mpr"][part]:
                    outputJson["region-wise-mpr"][part] = float('%.3g' % mprValue)
                if mprValue >= outputJsonMetric["value"]:
                    outputJsonMetric[part].append([float('%.3g' % float(mprLineSplit[1])), \
                    float('%.3g' % float(mprLineSplit[2])), float('%.3g' % float(mprLineSplit[3]))])
    jsonOutputFileMPr95 = 'MPr-95.json'
    with open(jsonOutputFileMPr95, 'w') as outfile:
        json.dump(outputJsonMetric, outfile, indent = 2, sort_keys=False)

# print(outputJson)
# jstr = json.dumps(outputJson, indent = 2, sort_keys=False)
# print(jstr)
# jsonOutputFile = uid + '_output_modified.json'
with open(jsonOutputFile, 'w') as outfile:
    json.dump(outputJson, outfile, indent = 2, sort_keys=False)
