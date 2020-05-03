# trace generated using paraview version 5.7.0-RC2
#
# To ensure correct image size when batch processing, please search 
# for and uncomment the line `# renderView*.ViewSize = [*,*]`

#### import the simple module from the paraview
from paraview.simple import *
import sys

uid = sys.argv[1]
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# get animation scene
animationScene1 = GetAnimationScene()

# get the time-keeper
timeKeeper1 = GetTimeKeeper()

# create a new 'PVD Reader'
coarse_brainpvd = PVDReader(FileName="coarse_brain_"+uid+".pvd")
coarse_brainpvd.CellArrays = ['PartID', 'AvgStrain', 'ProcID']
coarse_brainpvd.PointArrays = ['Displacements', 'Accelerations', 'Boundary']

# update animation scene based on data timesteps
animationScene1.UpdateAnimationUsingDataTimeSteps()

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')
# uncomment following to set a specific view size
# renderView1.ViewSize = [2174, 1719]

# show data in view
coarse_brainpvdDisplay = Show(coarse_brainpvd, renderView1)

# get color transfer function/color map for 'PartID'
partIDLUT = GetColorTransferFunction('PartID')

# get opacity transfer function/opacity map for 'PartID'
partIDPWF = GetOpacityTransferFunction('PartID')

# trace defaults for the display properties.
coarse_brainpvdDisplay.Representation = 'Surface'
coarse_brainpvdDisplay.ColorArrayName = ['CELLS', 'PartID']
coarse_brainpvdDisplay.LookupTable = partIDLUT
coarse_brainpvdDisplay.OSPRayScaleArray = 'Accelerations'
coarse_brainpvdDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
coarse_brainpvdDisplay.SelectOrientationVectors = 'Accelerations'
coarse_brainpvdDisplay.ScaleFactor = 0.0156247
coarse_brainpvdDisplay.SelectScaleArray = 'PartID'
coarse_brainpvdDisplay.GlyphType = 'Arrow'
coarse_brainpvdDisplay.GlyphTableIndexArray = 'PartID'
coarse_brainpvdDisplay.GaussianRadius = 0.000781235
coarse_brainpvdDisplay.SetScaleArray = ['POINTS', 'Accelerations']
coarse_brainpvdDisplay.ScaleTransferFunction = 'PiecewiseFunction'
coarse_brainpvdDisplay.OpacityArray = ['POINTS', 'Accelerations']
coarse_brainpvdDisplay.OpacityTransferFunction = 'PiecewiseFunction'
coarse_brainpvdDisplay.DataAxesGrid = 'GridAxesRepresentation'
coarse_brainpvdDisplay.PolarAxes = 'PolarAxesRepresentation'
coarse_brainpvdDisplay.ScalarOpacityFunction = partIDPWF
coarse_brainpvdDisplay.ScalarOpacityUnitDistance = 0.009666865901312131

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
coarse_brainpvdDisplay.ScaleTransferFunction.Points = [-55.8035665, 0.0, 0.5, 0.0, 59.251199, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
coarse_brainpvdDisplay.OpacityTransferFunction.Points = [-55.8035665, 0.0, 0.5, 0.0, 59.251199, 1.0, 0.5, 0.0]

# reset view to fit data
renderView1.ResetCamera()

# get the material library
materialLibrary1 = GetMaterialLibrary()

# show color bar/color legend
coarse_brainpvdDisplay.SetScalarBarVisibility(renderView1, True)

# update the view to ensure updated data information
renderView1.Update()

# create a new 'Clip'
clip1 = Clip(Input=coarse_brainpvd)
clip1.ClipType = 'Plane'
clip1.Scalars = ['CELLS', 'PartID']
clip1.Value = 4.5

# init the 'Plane' selected for 'ClipType'
clip1.ClipType.Origin = [0.11585740805, -0.3572862925, -0.0275435]

# toggle 3D widget visibility (only when running from the GUI)
Hide3DWidgets(proxy=clip1.ClipType)

# Properties modified on clip1
clip1.Invert = 0

# Properties modified on clip1.ClipType
clip1.ClipType.Normal = [0.0, 0.0, 1.0]

# show data in view
clip1Display = Show(clip1, renderView1)

# trace defaults for the display properties.
clip1Display.Representation = 'Surface'
clip1Display.ColorArrayName = ['CELLS', 'PartID']
clip1Display.LookupTable = partIDLUT
clip1Display.OSPRayScaleArray = 'Accelerations'
clip1Display.OSPRayScaleFunction = 'PiecewiseFunction'
clip1Display.SelectOrientationVectors = 'Accelerations'
clip1Display.ScaleFactor = 0.013361044099999997
clip1Display.SelectScaleArray = 'PartID'
clip1Display.GlyphType = 'Arrow'
clip1Display.GlyphTableIndexArray = 'PartID'
clip1Display.GaussianRadius = 0.0006680522049999999
clip1Display.SetScaleArray = ['POINTS', 'Accelerations']
clip1Display.ScaleTransferFunction = 'PiecewiseFunction'
clip1Display.OpacityArray = ['POINTS', 'Accelerations']
clip1Display.OpacityTransferFunction = 'PiecewiseFunction'
clip1Display.DataAxesGrid = 'GridAxesRepresentation'
clip1Display.PolarAxes = 'PolarAxesRepresentation'
clip1Display.ScalarOpacityFunction = partIDPWF
clip1Display.ScalarOpacityUnitDistance = 0.009401949892185164

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
clip1Display.ScaleTransferFunction.Points = [-55.8035665, 0.0, 0.5, 0.0, 59.251199, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
clip1Display.OpacityTransferFunction.Points = [-55.8035665, 0.0, 0.5, 0.0, 59.251199, 1.0, 0.5, 0.0]

# hide data in view
Hide(coarse_brainpvd, renderView1)

# show color bar/color legend
clip1Display.SetScalarBarVisibility(renderView1, True)

# update the view to ensure updated data information
renderView1.Update()

# set scalar coloring
ColorBy(clip1Display, ('CELLS', 'AvgStrain', 'Magnitude'))

# Hide the scalar bar for this color map if no visible data is colored by it.
HideScalarBarIfNotNeeded(partIDLUT, renderView1)

# rescale color and/or opacity maps used to include current data range
clip1Display.RescaleTransferFunctionToDataRange(True, False)

# show color bar/color legend
clip1Display.SetScalarBarVisibility(renderView1, True)

# get color transfer function/color map for 'AvgStrain'
avgStrainLUT = GetColorTransferFunction('AvgStrain')

# get opacity transfer function/opacity map for 'AvgStrain'
avgStrainPWF = GetOpacityTransferFunction('AvgStrain')

# current camera placement for renderView1
renderView1.CameraPosition = [0.1552542364739311, -0.38254250195511263, -0.49289536586670646]
renderView1.CameraFocalPoint = [0.1158574080500003, -0.35728629249999966, -0.027543499999999835]
renderView1.CameraViewUp = [0.5635994966584945, -0.820879116725902, 0.09226637030681739]
renderView1.CameraParallelScale = 0.1210494059939952

# save animation
SaveAnimation('simulation_'+uid+'.png', renderView1, ImageResolution=[1084, 856], FrameWindow=[0, 50])
