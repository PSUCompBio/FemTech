# generated using paraview version 5.8.0-636-gfeb3424073
#
# To ensure correct image size when batch processing, please search 
# for and uncomment the line `# renderView*.ViewSize = [*,*]`

#### import the simple module from the paraview
from paraview.simple import *
import sys

uid = sys.argv[1]
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

LoadPalette(paletteName='BlackBackground')

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')
# uncomment following to set a specific view size
# renderView1.ViewSize = [1154, 839]

# get layout
layout1 = GetLayout()

# get layout
layout1_1 = GetLayout()

# split cell
layout1_1.SplitHorizontal(0, 0.5)

# set active view
SetActiveView(None)

# Create a new 'Render View'
renderView2 = CreateView('RenderView')
renderView2.AxesGrid = 'GridAxes3DActor'
renderView2.StereoType = 'Crystal Eyes'
renderView2.CameraFocalDisk = 1.0
# uncomment following to set a specific view size
# renderView2.ViewSize = [400, 400]

# assign view to a particular cell in the layout
AssignViewToLayout(view=renderView2, layout=layout1_1, hint=2)

# split cell
layout1_1.SplitVertical(2, 0.5)

# set active view
SetActiveView(None)

# Create a new 'Render View'
renderView3 = CreateView('RenderView')
renderView3.AxesGrid = 'GridAxes3DActor'
renderView3.StereoType = 'Crystal Eyes'
renderView3.CameraFocalDisk = 1.0
# uncomment following to set a specific view size
# renderView3.ViewSize = [400, 400]

# assign view to a particular cell in the layout
AssignViewToLayout(view=renderView3, layout=layout1_1, hint=6)

# set active view
SetActiveView(renderView1)

# split cell
layout1_1.SplitVertical(1, 0.5)

# set active view
SetActiveView(None)

# Create a new 'Render View'
renderView4 = CreateView('RenderView')
renderView4.AxesGrid = 'GridAxes3DActor'
renderView4.StereoType = 'Crystal Eyes'
renderView4.CameraFocalDisk = 1.0
# uncomment following to set a specific view size
# renderView4.ViewSize = [400, 400]

# assign view to a particular cell in the layout
AssignViewToLayout(view=renderView4, layout=layout1_1, hint=4)

# get animation scene
animationScene1 = GetAnimationScene()

# get the time-keeper
timeKeeper1 = GetTimeKeeper()

# create a new 'PVD Reader'
coarse_brain_harry_s0JiuaYffpvd = PVDReader(FileName=uid+".pvd")
coarse_brain_harry_s0JiuaYffpvd.CellArrays = ['PartID', 'AvgStrain', 'ProcID']
coarse_brain_harry_s0JiuaYffpvd.PointArrays = ['Displacements', 'Accelerations', 'Boundary']

# update animation scene based on data timesteps
animationScene1.UpdateAnimationUsingDataTimeSteps()

# show data in view
coarse_brain_harry_s0JiuaYffpvdDisplay = Show(coarse_brain_harry_s0JiuaYffpvd, renderView4, 'UnstructuredGridRepresentation')

# get color transfer function/color map for 'PartID'
partIDLUT = GetColorTransferFunction('PartID')

# get opacity transfer function/opacity map for 'PartID'
partIDPWF = GetOpacityTransferFunction('PartID')

# trace defaults for the display properties.
coarse_brain_harry_s0JiuaYffpvdDisplay.Representation = 'Surface'
coarse_brain_harry_s0JiuaYffpvdDisplay.ColorArrayName = ['CELLS', 'PartID']
coarse_brain_harry_s0JiuaYffpvdDisplay.LookupTable = partIDLUT
coarse_brain_harry_s0JiuaYffpvdDisplay.OSPRayScaleArray = 'PartID'
coarse_brain_harry_s0JiuaYffpvdDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
coarse_brain_harry_s0JiuaYffpvdDisplay.SelectOrientationVectors = 'Accelerations'
coarse_brain_harry_s0JiuaYffpvdDisplay.ScaleFactor = 0.0156247
coarse_brain_harry_s0JiuaYffpvdDisplay.SelectScaleArray = 'PartID'
coarse_brain_harry_s0JiuaYffpvdDisplay.GlyphType = 'Arrow'
coarse_brain_harry_s0JiuaYffpvdDisplay.GlyphTableIndexArray = 'PartID'
coarse_brain_harry_s0JiuaYffpvdDisplay.GaussianRadius = 0.000781235
coarse_brain_harry_s0JiuaYffpvdDisplay.SetScaleArray = ['POINTS', 'PartID']
coarse_brain_harry_s0JiuaYffpvdDisplay.ScaleTransferFunction = 'PiecewiseFunction'
coarse_brain_harry_s0JiuaYffpvdDisplay.OpacityArray = ['POINTS', 'PartID']
coarse_brain_harry_s0JiuaYffpvdDisplay.OpacityTransferFunction = 'PiecewiseFunction'
coarse_brain_harry_s0JiuaYffpvdDisplay.DataAxesGrid = 'GridAxesRepresentation'
coarse_brain_harry_s0JiuaYffpvdDisplay.PolarAxes = 'PolarAxesRepresentation'
coarse_brain_harry_s0JiuaYffpvdDisplay.ScalarOpacityFunction = partIDPWF
coarse_brain_harry_s0JiuaYffpvdDisplay.ScalarOpacityUnitDistance = 0.009188415966156702

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
coarse_brain_harry_s0JiuaYffpvdDisplay.ScaleTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 9.0, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
coarse_brain_harry_s0JiuaYffpvdDisplay.OpacityTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 9.0, 1.0, 0.5, 0.0]

# reset view to fit data
renderView4.ResetCamera()

# show color bar/color legend
coarse_brain_harry_s0JiuaYffpvdDisplay.SetScalarBarVisibility(renderView4, True)

# update the view to ensure updated data information
renderView4.Update()

# create a new 'Threshold'
threshold1 = Threshold(Input=coarse_brain_harry_s0JiuaYffpvd)
threshold1.Scalars = ['POINTS', 'PartID']
threshold1.ThresholdRange = [0.0, 9.0]

# show data in view
threshold1Display = Show(threshold1, renderView4, 'UnstructuredGridRepresentation')

# trace defaults for the display properties.
threshold1Display.Representation = 'Surface'
threshold1Display.ColorArrayName = ['CELLS', 'PartID']
threshold1Display.LookupTable = partIDLUT
threshold1Display.OSPRayScaleArray = 'PartID'
threshold1Display.OSPRayScaleFunction = 'PiecewiseFunction'
threshold1Display.SelectOrientationVectors = 'Accelerations'
threshold1Display.ScaleFactor = 0.0156247
threshold1Display.SelectScaleArray = 'PartID'
threshold1Display.GlyphType = 'Arrow'
threshold1Display.GlyphTableIndexArray = 'PartID'
threshold1Display.GaussianRadius = 0.000781235
threshold1Display.SetScaleArray = ['POINTS', 'PartID']
threshold1Display.ScaleTransferFunction = 'PiecewiseFunction'
threshold1Display.OpacityArray = ['POINTS', 'PartID']
threshold1Display.OpacityTransferFunction = 'PiecewiseFunction'
threshold1Display.DataAxesGrid = 'GridAxesRepresentation'
threshold1Display.PolarAxes = 'PolarAxesRepresentation'
threshold1Display.ScalarOpacityFunction = partIDPWF
threshold1Display.ScalarOpacityUnitDistance = 0.009188415966156702

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
threshold1Display.ScaleTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 9.0, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
threshold1Display.OpacityTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 9.0, 1.0, 0.5, 0.0]

# hide data in view
Hide(coarse_brain_harry_s0JiuaYffpvd, renderView4)

# show color bar/color legend
threshold1Display.SetScalarBarVisibility(renderView4, True)

# update the view to ensure updated data information
renderView4.Update()

# Properties modified on threshold1
threshold1.ThresholdRange = [1.0, 9.0]

# update the view to ensure updated data information
renderView4.Update()

# change representation type
threshold1Display.SetRepresentationType('Volume')

# set scalar coloring
ColorBy(threshold1Display, ('POINTS', 'Boundary', 'Magnitude'))

# Hide the scalar bar for this color map if no visible data is colored by it.
HideScalarBarIfNotNeeded(partIDLUT, renderView4)

# rescale color and/or opacity maps used to include current data range
threshold1Display.RescaleTransferFunctionToDataRange(True, False)

# show color bar/color legend
threshold1Display.SetScalarBarVisibility(renderView4, True)

# get color transfer function/color map for 'Boundary'
boundaryLUT = GetColorTransferFunction('Boundary')

# get opacity transfer function/opacity map for 'Boundary'
boundaryPWF = GetOpacityTransferFunction('Boundary')

# hide color bar/color legend
threshold1Display.SetScalarBarVisibility(renderView4, False)

# Apply a preset using its name. Note this may not work as expected when presets have duplicate names.
boundaryLUT.ApplyPreset('Blue Orange (divergent)', True)

# Properties modified on threshold1
threshold1.AllScalars = 0

# update the view to ensure updated data information
renderView4.Update()

# reset view to fit data
renderView4.ResetCamera()

# Add new text box
annotateTimeFilter2 = AnnotateTimeFilter(Input=threshold1)
annotateTimeFilter2Display = Show(annotateTimeFilter2, renderView1, 'TextSourceRepresentation')
annotateTimeFilter2.Format = 'Axes : J211'
annotateTimeFilter2Display.WindowLocation = 'AnyLocation'
annotateTimeFilter2Display.Position = [0.01, 0.8]
renderView1.Update()

# create a new 'Annotate Time Filter'
annotateTimeFilter1 = AnnotateTimeFilter(Input=threshold1)

# show data in view
annotateTimeFilter1Display = Show(annotateTimeFilter1, renderView4, 'TextSourceRepresentation')

# update the view to ensure updated data information
renderView4.Update()

# Properties modified on annotateTimeFilter1
annotateTimeFilter1.Format = 'Time: %6.2f (msec)'

# update the view to ensure updated data information
renderView4.Update()

# Properties modified on annotateTimeFilter1
annotateTimeFilter1.Scale = 1000.0

# update the view to ensure updated data information
renderView4.Update()

# hide data in view
Hide(annotateTimeFilter1, renderView4)

# set active view
SetActiveView(renderView1)

# set active source
SetActiveSource(annotateTimeFilter1)

# show data in view
annotateTimeFilter1Display_1 = Show(annotateTimeFilter1, renderView1, 'TextSourceRepresentation')

# reset view to fit data
renderView1.ResetCamera()

# set active source
SetActiveSource(threshold1)

# show data in view
threshold1Display_1 = Show(threshold1, renderView1, 'UnstructuredGridRepresentation')

# trace defaults for the display properties.
threshold1Display_1.Representation = 'Surface'
threshold1Display_1.ColorArrayName = ['CELLS', 'PartID']
threshold1Display_1.LookupTable = partIDLUT
threshold1Display_1.OSPRayScaleArray = 'PartID'
threshold1Display_1.OSPRayScaleFunction = 'PiecewiseFunction'
threshold1Display_1.SelectOrientationVectors = 'Accelerations'
threshold1Display_1.ScaleFactor = 0.015142000000000001
threshold1Display_1.SelectScaleArray = 'PartID'
threshold1Display_1.GlyphType = 'Arrow'
threshold1Display_1.GlyphTableIndexArray = 'PartID'
threshold1Display_1.GaussianRadius = 0.0007571
threshold1Display_1.SetScaleArray = ['POINTS', 'PartID']
threshold1Display_1.ScaleTransferFunction = 'PiecewiseFunction'
threshold1Display_1.OpacityArray = ['POINTS', 'PartID']
threshold1Display_1.OpacityTransferFunction = 'PiecewiseFunction'
threshold1Display_1.DataAxesGrid = 'GridAxesRepresentation'
threshold1Display_1.PolarAxes = 'PolarAxesRepresentation'
threshold1Display_1.ScalarOpacityFunction = partIDPWF
threshold1Display_1.ScalarOpacityUnitDistance = 0.009561245734043766

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
threshold1Display_1.ScaleTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 9.0, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
threshold1Display_1.OpacityTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 9.0, 1.0, 0.5, 0.0]

# show color bar/color legend
threshold1Display_1.SetScalarBarVisibility(renderView1, True)

# update the view to ensure updated data information
renderView1.Update()

# reset view to fit data bounds
renderView1.ResetCamera(-0.082497, 0.068923, -0.057284, 0.057284, -0.088786, 0.031734)

# reset view to fit data
renderView1.ResetCamera()

# change representation type
threshold1Display_1.SetRepresentationType('Volume')

# set scalar coloring
ColorBy(threshold1Display_1, ('POINTS', 'Boundary', 'Magnitude'))

# Hide the scalar bar for this color map if no visible data is colored by it.
HideScalarBarIfNotNeeded(partIDLUT, renderView1)

# rescale color and/or opacity maps used to include current data range
threshold1Display_1.RescaleTransferFunctionToDataRange(True, False)

# show color bar/color legend
threshold1Display_1.SetScalarBarVisibility(renderView1, True)

# hide color bar/color legend
threshold1Display_1.SetScalarBarVisibility(renderView1, False)

# set active view
SetActiveView(renderView2)

# show data in view
threshold1Display_2 = Show(threshold1, renderView2, 'UnstructuredGridRepresentation')

# trace defaults for the display properties.
threshold1Display_2.Representation = 'Surface'
threshold1Display_2.ColorArrayName = ['CELLS', 'PartID']
threshold1Display_2.LookupTable = partIDLUT
threshold1Display_2.OSPRayScaleArray = 'PartID'
threshold1Display_2.OSPRayScaleFunction = 'PiecewiseFunction'
threshold1Display_2.SelectOrientationVectors = 'Accelerations'
threshold1Display_2.ScaleFactor = 0.015142000000000001
threshold1Display_2.SelectScaleArray = 'PartID'
threshold1Display_2.GlyphType = 'Arrow'
threshold1Display_2.GlyphTableIndexArray = 'PartID'
threshold1Display_2.GaussianRadius = 0.0007571
threshold1Display_2.SetScaleArray = ['POINTS', 'PartID']
threshold1Display_2.ScaleTransferFunction = 'PiecewiseFunction'
threshold1Display_2.OpacityArray = ['POINTS', 'PartID']
threshold1Display_2.OpacityTransferFunction = 'PiecewiseFunction'
threshold1Display_2.DataAxesGrid = 'GridAxesRepresentation'
threshold1Display_2.PolarAxes = 'PolarAxesRepresentation'
threshold1Display_2.ScalarOpacityFunction = partIDPWF
threshold1Display_2.ScalarOpacityUnitDistance = 0.009561245734043766

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
threshold1Display_2.ScaleTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 9.0, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
threshold1Display_2.OpacityTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 9.0, 1.0, 0.5, 0.0]

# show color bar/color legend
threshold1Display_2.SetScalarBarVisibility(renderView2, True)

# reset view to fit data
renderView2.ResetCamera()

# reset view to fit data
renderView2.ResetCamera()

# change representation type
threshold1Display_2.SetRepresentationType('Volume')

# set scalar coloring
ColorBy(threshold1Display_2, ('POINTS', 'Boundary', 'Magnitude'))

# Hide the scalar bar for this color map if no visible data is colored by it.
HideScalarBarIfNotNeeded(partIDLUT, renderView2)

# rescale color and/or opacity maps used to include current data range
threshold1Display_2.RescaleTransferFunctionToDataRange(True, False)

# show color bar/color legend
threshold1Display_2.SetScalarBarVisibility(renderView2, True)

# hide color bar/color legend
threshold1Display_2.SetScalarBarVisibility(renderView2, False)

# set active view
SetActiveView(renderView1)

# set active view
SetActiveView(renderView4)

# set active view
SetActiveView(renderView2)

# current camera placement for renderView1
renderView1.CameraPosition = [0.9007775473613021, 0.0, -0.028526000000000003]
renderView1.CameraFocalPoint = [0.46630882042000044, 0.0, -0.028526000000000003]
renderView1.CameraViewUp = [0.0, 4.440892098500626e-16, -1.0]
renderView1.CameraParallelScale = 0.11244878103385558

# current camera placement for renderView2
renderView2.CameraPosition = [-0.006787000000000001, 0.0, -1.1683530781692515]
renderView2.CameraFocalPoint = [-0.006787000000000001, 0.0, -0.73388435122795]
renderView2.CameraViewUp = [-1.0, 2.220446049250313e-16, 0.0]
renderView2.CameraParallelScale = 0.11244878103385558

# current camera placement for renderView3

# current camera placement for renderView4
renderView4.CameraPosition = [-0.006787000000000001, -0.9979077462412631, -0.028526000000000003]
renderView4.CameraFocalPoint = [-0.006787000000000001, -0.5634390192999611, -0.028526000000000003]
renderView4.CameraViewUp = [4.440892098500626e-16, 0.0, -1.0]
renderView4.CameraParallelScale = 0.11244878103385558

# set active view
SetActiveView(renderView3)
renderView3.OrientationAxesVisibility = 0

# save animation
SaveAnimation('simulation_'+uid+'.png', layout1_1, SaveAllViews=1,
    ImageResolution=[1920, 1080], FrameWindow=[0, 50])
