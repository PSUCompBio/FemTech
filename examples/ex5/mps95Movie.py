# trace generated using paraview version 5.8.0
#
# To ensure correct image size when batch processing, please search 
# for and uncomment the line `# renderView*.ViewSize = [*,*]`

#### import the simple module from the paraview
from paraview.simple import *
import sys
import json
import os
import numpy as np

jsonFile = sys.argv[1]
inputJson = json.loads(open(jsonFile).read())
uid = inputJson["uid"]
simulationJson = inputJson["simulation"]
meshFile = simulationJson["mesh"]
mesh = os.path.splitext(meshFile)[0]

# Read time steps reported in pvd files
paraviewFile = mesh+"_"+uid+".pvd"

#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

defaultCameraPosition = [-0.2007737695527859, -0.16394057453798055, -0.5569876615568349]
defaultCameraFocalPoint = [-0.05191673062823166, -0.06371209375912919, -0.13701646316140037]
defaultCameraViewUp = [0.05342884714081945, -0.9754046954694863, 0.21384816657919195]

# defaultCameraPosition = [0.05887805261889896, -0.4059814608821956, 0.1278711190042479]
# defaultCameraFocalPoint = [0.02730430718751414, 0.029470759431911277, -0.006160356035267402]
# defaultCameraViewUp = [0.8853796124608853, 0.19419934515902124, 0.4223618782257934]
defaultCameraParallelScale = 0.11820409844628908

offset = np.zeros(3)
headCG = np.zeros(3)
if "angular-sensor-position" in inputJson["simulation"]:
    offset = np.array(inputJson["simulation"]["angular-sensor-position"])
    # get CG
    if "head-cg" in inputJson["simulation"]:
        headCG = np.array(inputJson["simulation"]["head-cg"])
    # else:
        # Get CG from log
offset = offset-headCG
defaultCameraFocalPoint = defaultCameraFocalPoint - offset
defaultCameraPosition = defaultCameraPosition - offset

if "mesh-transformation" in inputJson["simulation"]:
    transform = inputJson["simulation"]["mesh-transformation"]
    signP = np.zeros(3)
    mapP = np.zeros(3, dtype=int)
    for i, st in enumerate(transform):
        if len(st) == 1:
            signP[i] = 1;
            if st[0] == 'x':
                mapP[i] = 0
            else:
                if st[0] == 'y':
                    mapP[i] = 1
                else:
                    if st[0] == 'z':
                        mapP[i] = 2
                    else:
                        print('Error in mesh transformation')
        else:
            if len(st) == 2:
                if st[0] == '-':
                    signP[i] = -1
                else:
                    if st[0] == '+':
                        signP[i] = 1
                    else:
                        print('Error in mesh transformation array sign')
                if st[1] == 'x':
                    mapP[i] = 0
                else:
                    if st[1] == 'y':
                        mapP[i] = 1
                    else:
                        if st[1] == 'z':
                            mapP[i] = 2
                        else:
                            print('Error in mesh transformation')
            else:
                print('Error in mesh transformation array string')

    transformedCameraPosition = np.zeros(3)
    transformedCameraFocalPoint = np.zeros(3)
    transformedCameraViewUp = np.zeros(3)

    for i in np.arange(3):
        transformedCameraPosition[mapP[i]] = defaultCameraPosition[i]*signP[i]
        transformedCameraFocalPoint[mapP[i]] = defaultCameraFocalPoint[i]*signP[i]
        transformedCameraViewUp[mapP[i]] = defaultCameraViewUp[i]*signP[i]
else:
    transformedCameraPosition = defaultCameraPosition
    transformedCameraFocalPoint = defaultCameraFocalPoint
    transformedCameraViewUp = defaultCameraViewUp

# print(transformedCameraPosition)
# print(transformedCameraFocalPoint)
# print(transformedCameraViewUp)
#
# create a new 'PVD Reader'
coarse_brain_test_ptpvd = PVDReader(FileName=paraviewFile)
coarse_brain_test_ptpvd.CellArrays = ['PartID', 'AvgStrain', 'ProcID', 'CSDM-5', 'CSDM-10', 'CSDM-15', 'CSDM-30', 'MPSR-120', 'MPSxSR-28', 'MPS-95']
coarse_brain_test_ptpvd.PointArrays = ['Displacements', 'Accelerations', 'Boundary']

# get animation scene
animationScene1 = GetAnimationScene()
# get active view
renderView1 = GetActiveViewOrCreate('RenderView')
# uncomment following to set a specific view size
renderView1.ViewSize = [1538, 838]

# get color transfer function/color map for 'PartID'
partIDLUT = GetColorTransferFunction('PartID')

# update animation scene based on data timesteps
animationScene1.UpdateAnimationUsingDataTimeSteps()

# create a new 'Threshold'
threshold1 = Threshold(Input=coarse_brain_test_ptpvd)
threshold1.Scalars = ['POINTS', 'PartID']
threshold1.ThresholdRange = [2.0, 9.0]

# show data in view
threshold1Display = Show(threshold1, renderView1, 'UnstructuredGridRepresentation')

Hide(coarse_brain_test_ptpvd, renderView1)

# show color bar/color legend
threshold1Display.SetScalarBarVisibility(renderView1, True)
# hide color bar/color legend
threshold1Display.SetScalarBarVisibility(renderView1, False)

# turn off scalar coloring
ColorBy(threshold1Display, None)

# Hide the scalar bar for this color map if no visible data is colored by it.
HideScalarBarIfNotNeeded(partIDLUT, renderView1)

# change solid color
threshold1Display.AmbientColor = [0.0, 0.6666666666666666, 1.0]
threshold1Display.DiffuseColor = [0.0, 0.6666666666666666, 1.0]

# Properties modified on threshold1Display
threshold1Display.Opacity = 0.1

# create a new 'Threshold'
threshold2 = Threshold(Input=threshold1)
threshold2.Scalars = ['POINTS', 'MPS-95']
threshold2.ThresholdRange = [1.0, 1.0]
threshold2.AllScalars = 0

# show data in view
threshold2Display = Show(threshold2, renderView1, 'UnstructuredGridRepresentation')

# trace defaults for the display properties.
threshold2Display.Representation = 'Surface'
threshold2Display.ColorArrayName = ['CELLS', 'PartID']
# set scalar coloring
ColorBy(threshold2Display, ('CELLS', 'PartID'))

# show color bar/color legend
threshold2Display.SetScalarBarVisibility(renderView1, True)

# Apply a preset using its name. Note this may not work as expected when presets have duplicate names.
partIDLUT.ApplyPreset('Viridis (matplotlib)', True)

# get color legend/bar for partIDLUT in view renderView1
partIDLUTColorBar = GetScalarBar(partIDLUT, renderView1)

partIDLUTColorBar.WindowLocation = 'AnyLocation'
partIDLUTColorBar.Position = [0.9057217165149545, 0.1966626936829559]
partIDLUTColorBar.ScalarBarLength = 0.3299999999999999
partIDLUTColorBar.ScalarBarThickness = 16
partIDLUTColorBar.TitleFontFamily = 'Times'
partIDLUTColorBar.TitleFontSize = 14
partIDLUTColorBar.LabelFontFamily = 'Times'
partIDLUTColorBar.LabelFontSize = 14
partIDLUTColorBar.AutomaticLabelFormat = 1
partIDLUTColorBar.AddRangeLabels = 0

annotateTime1 = AnnotateTime()
annotateTimeFilter1 = AnnotateTimeFilter(Input=annotateTime1)
annotateTimeFilter1Display = Show(annotateTimeFilter1, renderView1, 'TextSourceRepresentation')
annotateTimeFilter1.Format = 'Time: %6.2f (msec)'
annotateTimeFilter1.Scale = 1000.0
annotateTimeFilter1Display.FontFamily = 'Times'
annotateTimeFilter1Display.WindowLocation = 'AnyLocation'
annotateTimeFilter1Display.FontSize = 10
annotateTimeFilter1Display.Position = [0.01, 0.94]

# current camera placement for renderView1
renderView1.CameraPosition = transformedCameraPosition
renderView1.CameraFocalPoint = transformedCameraFocalPoint
renderView1.CameraViewUp = transformedCameraViewUp
renderView1.CameraParallelScale = defaultCameraParallelScale

# save animation
SaveAnimation('injury_'+uid+'.png', renderView1, ImageResolution=[1538, 838])
