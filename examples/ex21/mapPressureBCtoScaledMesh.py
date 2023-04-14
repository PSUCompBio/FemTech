import json

import numpy as np

pressureInput = "pressureBCRaw.json"
pressureMapFile = "pressureToNodeIDMap.json"

pressureJson = json.loads(open(pressureInput).read())
pressureMap = json.loads(open(pressureMapFile).read())

nPressure = len(pressureJson)
nNodes = len(pressureMap)

nodeID = np.zeros(nNodes, dtype=int)
pressureVec = np.zeros(nNodes)

for i in np.arange(0, nNodes):
    node = pressureMap[i]

for i in np.arange(0, nNodes):
    elem = pressureMap[i]
    nodeID[i] = int(elem['nodeID'])
    pressureID = int(elem['pressureID'])
    pressureVec[i] = float(pressureJson[pressureID]['pressure_value'])

print(nodeID)
print(pressureVec)

jsonObject = []
for i in np.arange(0, nNodes):
    elem = {}
    elem['nodeID'] = int(nodeID[i])
    elem['pressure'] = pressureVec[i]
    jsonObject.append(elem)

with open('pressureBC.json', 'w') as f:
    json.dump(jsonObject, f)

