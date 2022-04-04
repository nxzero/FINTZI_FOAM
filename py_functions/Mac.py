# trace generated using paraview version 5.9.1

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# get active source.
a0025foam = GetActiveSource()

# Properties modified on a0025foam
a0025foam.MeshRegions = ['air_to_sphere', 'air_to_spherePlaneFacesbot', 'air_to_spherePlaneFacestop']

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')

# show data in view
a0025foamDisplay = Show(a0025foam, renderView1, 'GeometryRepresentation')

# get color transfer function/color map for 'p'
pLUT = GetColorTransferFunction('p')

# trace defaults for the display properties.
a0025foamDisplay.Representation = 'Surface'
a0025foamDisplay.ColorArrayName = ['POINTS', 'p']
a0025foamDisplay.LookupTable = pLUT
a0025foamDisplay.SelectTCoordArray = 'None'
a0025foamDisplay.SelectNormalArray = 'None'
a0025foamDisplay.SelectTangentArray = 'None'
a0025foamDisplay.OSPRayScaleArray = 'p'
a0025foamDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
a0025foamDisplay.SelectOrientationVectors = 'U'
a0025foamDisplay.ScaleFactor = 0.166064453125
a0025foamDisplay.SelectScaleArray = 'p'
a0025foamDisplay.GlyphType = 'Arrow'
a0025foamDisplay.GlyphTableIndexArray = 'p'
a0025foamDisplay.GaussianRadius = 0.00830322265625
a0025foamDisplay.SetScaleArray = ['POINTS', 'p']
a0025foamDisplay.ScaleTransferFunction = 'PiecewiseFunction'
a0025foamDisplay.OpacityArray = ['POINTS', 'p']
a0025foamDisplay.OpacityTransferFunction = 'PiecewiseFunction'
a0025foamDisplay.DataAxesGrid = 'GridAxesRepresentation'
a0025foamDisplay.PolarAxes = 'PolarAxesRepresentation'

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
a0025foamDisplay.ScaleTransferFunction.Points = [-1.4575281143188477, 0.0, 0.5, 0.0, 0.8620374798774719, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
a0025foamDisplay.OpacityTransferFunction.Points = [-1.4575281143188477, 0.0, 0.5, 0.0, 0.8620374798774719, 1.0, 0.5, 0.0]

# reset view to fit data
renderView1.ResetCamera()

# show color bar/color legend
a0025foamDisplay.SetScalarBarVisibility(renderView1, True)

# update the view to ensure updated data information
renderView1.Update()

# get opacity transfer function/opacity map for 'p'
pPWF = GetOpacityTransferFunction('p')

# create a new 'Integrate Variables'
integrateVariables1 = IntegrateVariables(registrationName='IntegrateVariables1', Input=a0025foam)

# Create a new 'SpreadSheet View'
spreadSheetView1 = CreateView('SpreadSheetView')
spreadSheetView1.ColumnToSort = ''
spreadSheetView1.BlockSize = 1024

# show data in view
integrateVariables1Display = Show(integrateVariables1, spreadSheetView1, 'SpreadSheetRepresentation')

# get layout
layout1 = GetLayoutByName("Layout #1")

# add view to a layout so it's visible in UI
AssignViewToLayout(view=spreadSheetView1, layout=layout1, hint=0)

# show data in view
a0025foamDisplay_1 = Show(a0025foam, spreadSheetView1, 'SpreadSheetRepresentation')

# get animation scene
animationScene1 = GetAnimationScene()

animationScene1.GoToLast()

# export view
ExportView('datas.csv', view=spreadSheetView1)

#================================================================
# addendum: following script captures some of the application
# state to faithfully reproduce the visualization during playback
#================================================================

#--------------------------------
# saving layout sizes for layouts

# layout/tab size in pixels
layout1.SetSize(1104, 781)

#-----------------------------------
# saving camera placements for views

# current camera placement for renderView1
renderView1.CameraPosition = [17.21122804909539, 18.32059076707465, 20.913440521075028]
renderView1.CameraFocalPoint = [17.499999999999996, 17.499999999999996, 17.499999999999996]
renderView1.CameraViewUp = [0.09684561543150856, 0.9695596141370174, -0.2248890424324192]
renderView1.CameraParallelScale = 1.1031596022115961

#--------------------------------------------
# uncomment the following to render all views
# RenderAllViews()
# alternatively, if you want to write images, you can use SaveScreenshot(...).