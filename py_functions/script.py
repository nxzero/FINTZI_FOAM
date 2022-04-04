# trace generated using paraview version 5.9.1

#### import the simple module from the paraview
from paraview.simple import *
import numpy as np
import math
import sys
sys.path.append('.')
import parameters_for_this_study as p
import os
D = p.parameters['r']*2
L = p.parameters['length']
chi = p.parameters['ksi']
xCenter =p.parameters['xCenter']
yCenter =p.parameters['yCenter']
zCenter =p.parameters['zCenter']
dH  = D/4. 
N = int(L/(D/4.))
ThetaU = p.parameters['ThetaU']
Theta = p.parameters['Theta']
ct = np.cos(ThetaU*math.pi/180)
ctt = np.cos(Theta*math.pi/180)
st = np.sin(ThetaU*math.pi/180)
stt = np.sin(Theta*math.pi/180)
Re = p.parameters['Re']
HOME = os.getenv('HOME')
Cond = str(ctt)+'*(coordsX-'+str(xCenter)+')   +   '+str(stt)+'*(coordsY-'+str(yCenter)+')'
nx = ctt
ny = stt

os.system('mkdir '+HOME+'/CYLINDERS_PROJECT')
os.system('mkdir '+HOME+'/CYLINDERS_PROJECT/pp')
os.system('mkdir '+HOME+'/CYLINDERS_PROJECT/pp/forceslin')
os.system('mkdir '+HOME+'/CYLINDERS_PROJECT/pp/forceslin/Chi'+str(chi))
os.system('mkdir '+HOME+'/CYLINDERS_PROJECT/pp/forceslin/Chi'+str(chi)+'/Re'+str(Re))
os.system('mkdir '+HOME+'/CYLINDERS_PROJECT/pp/forceslin/Chi'+str(chi)+'/Re'+str(Re)+'/FD')
os.system('mkdir '+HOME+'/CYLINDERS_PROJECT/pp/forceslin/Chi'+str(chi)+'/Re'+str(Re)+'/FL')
os.system('mkdir '+HOME+'/CYLINDERS_PROJECT/pp/forceslin/Chi'+str(chi)+'/Re'+str(Re)+'/M')
PATH = HOME+'/CYLINDERS_PROJECT/pp/forceslin/Chi'+str(chi)+'/Re'+str(Re)+'/'
#### disable automatic camera reset on 'Show'

paraview.simple._DisableFirstRenderCameraReset()

# get active source.
a0025foam = GetActiveSource()

# Properties modified on a0025foam
a0025foam.MeshRegions = ['air_to_sphere']
a0025foam.CellArrays = ['U', 'forces:force', 'forces:moment', 'p']

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
a0025foamDisplay.ScaleFactor = 0.10000629425048829
a0025foamDisplay.SelectScaleArray = 'p'
a0025foamDisplay.GlyphType = 'Arrow'
a0025foamDisplay.GlyphTableIndexArray = 'p'
a0025foamDisplay.GaussianRadius = 0.005000314712524414
a0025foamDisplay.SetScaleArray = ['POINTS', 'p']
a0025foamDisplay.ScaleTransferFunction = 'PiecewiseFunction'
a0025foamDisplay.OpacityArray = ['POINTS', 'p']
a0025foamDisplay.OpacityTransferFunction = 'PiecewiseFunction'
a0025foamDisplay.DataAxesGrid = 'GridAxesRepresentation'
a0025foamDisplay.PolarAxes = 'PolarAxesRepresentation'

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
a0025foamDisplay.ScaleTransferFunction.Points = [-105.09708404541016, 0.0, 0.5, 0.0, 103.13551330566406, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
a0025foamDisplay.OpacityTransferFunction.Points = [-105.09708404541016, 0.0, 0.5, 0.0, 103.13551330566406, 1.0, 0.5, 0.0]

# reset view to fit data
renderView1.ResetCamera()

# show color bar/color legend
a0025foamDisplay.SetScalarBarVisibility(renderView1, True)

# update the view to ensure updated data information
renderView1.Update()
animationScene1 = GetAnimationScene()

animationScene1.GoToLast()

# get opacity transfer function/opacity map for 'p'
pPWF = GetOpacityTransferFunction('p')

# set scalar coloring
ColorBy(a0025foamDisplay, ('POINTS', 'forces:force', 'Magnitude'))

# Hide the scalar bar for this color map if no visible data is colored by it.
HideScalarBarIfNotNeeded(pLUT, renderView1)

# rescale color and/or opacity maps used to include current data range
a0025foamDisplay.RescaleTransferFunctionToDataRange(True, False)

# show color bar/color legend
a0025foamDisplay.SetScalarBarVisibility(renderView1, True)

# get color transfer function/color map for 'forcesforce'
forcesforceLUT = GetColorTransferFunction('forcesforce')

# get opacity transfer function/opacity map for 'forcesforce'
forcesforcePWF = GetOpacityTransferFunction('forcesforce')







for c in list(np.linspace(-L/2+dH/2,L/2-dH/2,20)):#list(np.linspace(-0.5,0.5,11)):
    # create a new 'Calculator'
    calculator1 = Calculator(registrationName='Calculator1', Input=a0025foam)
    calculator1.Function = ''
    # calculator1.ResultArrayName = 'Resul'+str(c)

    # Properties modified on calculator1
    C = c
    print('(forces:force_X*'+str(ct)+'+'+str(st)+'*forces:force_Y)*('+Cond+'<'+str(C+dH/2.)+')*('+Cond+'>'+str(C-dH/2.)+')')
    calculator1.Function = '(forces:force_X*'+str(ct)+'+'+str(st)+'*forces:force_Y)*('+Cond+'<'+str(C+dH/2.)+')*('+Cond+'>'+str(C-dH/2.)+')'

    # show data in view
    calculator1Display = Show(calculator1, renderView1, 'GeometryRepresentation')

    # get color transfer function/color map for 'Result'
    resultLUT = GetColorTransferFunction('Result')

    # trace defaults for the display properties.
    calculator1Display.Representation = 'Surface'
    calculator1Display.ColorArrayName = ['POINTS', 'Result']
    calculator1Display.LookupTable = resultLUT
    calculator1Display.SelectTCoordArray = 'None'
    calculator1Display.SelectNormalArray = 'None'
    calculator1Display.SelectTangentArray = 'None'
    calculator1Display.OSPRayScaleArray = 'Result'
    calculator1Display.OSPRayScaleFunction = 'PiecewiseFunction'
    calculator1Display.SelectOrientationVectors = 'U'
    calculator1Display.ScaleFactor = 0.10000629425048829
    calculator1Display.SelectScaleArray = 'Result'
    calculator1Display.GlyphType = 'Arrow'
    calculator1Display.GlyphTableIndexArray = 'Result'
    calculator1Display.GaussianRadius = 0.005000314712524414
    calculator1Display.SetScaleArray = ['POINTS', 'Result']
    calculator1Display.ScaleTransferFunction = 'PiecewiseFunction'
    calculator1Display.OpacityArray = ['POINTS', 'Result']
    calculator1Display.OpacityTransferFunction = 'PiecewiseFunction'
    calculator1Display.DataAxesGrid = 'GridAxesRepresentation'
    calculator1Display.PolarAxes = 'PolarAxesRepresentation'

    # init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
    calculator1Display.ScaleTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 14.734161467233076, 1.0, 0.5, 0.0]

    # init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
    calculator1Display.OpacityTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 14.734161467233076, 1.0, 0.5, 0.0]

    # hide data in view
    Hide(a0025foam, renderView1)

    # show color bar/color legend
    calculator1Display.SetScalarBarVisibility(renderView1, True)

    # update the view to ensure updated data information
    renderView1.Update()
    animationScene1 = GetAnimationScene()

    animationScene1.GoToLast()

    # get opacity transfer function/opacity map for 'Result'
    resultPWF = GetOpacityTransferFunction('Result')

    # create a new 'Integrate Variables'
    integrateVariables1 = IntegrateVariables(registrationName='IntegrateVariables1', Input=calculator1)
    # Properties modified on integrateVariables1
    # integrateVariables1.DivideCellDataByVolume = 1

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

    # Properties modified on spreadSheetView1
    spreadSheetView1.ColumnToSort = 'Result'
    spreadSheetView1.InvertOrder = 1

    # Properties modified on spreadSheetView1
    spreadSheetView1.InvertOrder = 0

    # Properties modified on spreadSheetView1
    spreadSheetView1.InvertOrder = 1

    # export view
    ExportView(PATH+'FD/forces'+str(round(c,4)).replace('.','_')+'.csv', view=spreadSheetView1)
    
    # destroy spreadSheetView1
    Delete(spreadSheetView1)
    del spreadSheetView1

    # close an empty frame
    layout1.Collapse(2)

    # set active view
    SetActiveView(renderView1)

    # set active source
    SetActiveSource(calculator1)

    # destroy integrateVariables1
    Delete(integrateVariables1)
    del integrateVariables1

    # set active source
    SetActiveSource(a0025foam)

    # hide data in view
    Hide(calculator1, renderView1)

    # show data in view
    a0025foamDisplay = Show(a0025foam, renderView1, 'GeometryRepresentation')

    # show color bar/color legend
    a0025foamDisplay.SetScalarBarVisibility(renderView1, True)

    # destroy calculator1
    Delete(calculator1)
    del calculator1
    ########################
    ########################
    ########################
    ########################
    ########################
    ########################
    ########################
    ########################
    ########################
    ########################
    ########################
    ########################
    ########################





for c in list(np.linspace(-L/2+dH/2,L/2-dH/2,20)):#list(np.linspace(-0.5,0.5,11)):
    # create a new 'Calculator'
    calculator1 = Calculator(registrationName='Calculator1', Input=a0025foam)
    calculator1.Function = ''
    # calculator1.ResultArrayName = 'Resul'+str(c)

    # Properties modified on calculator1
    C = c
    calculator1.Function = '(-forces:force_X*'+str(st)+'+'+str(ct)+'*forces:force_Y)*('+Cond+'<'+str(C+dH/2.)+')*('+Cond+'>'+str(C-dH/2.)+')'

    # show data in view
    calculator1Display = Show(calculator1, renderView1, 'GeometryRepresentation')

    # get color transfer function/color map for 'Result'
    resultLUT = GetColorTransferFunction('Result')

    # trace defaults for the display properties.
    calculator1Display.Representation = 'Surface'
    calculator1Display.ColorArrayName = ['POINTS', 'Result']
    calculator1Display.LookupTable = resultLUT
    calculator1Display.SelectTCoordArray = 'None'
    calculator1Display.SelectNormalArray = 'None'
    calculator1Display.SelectTangentArray = 'None'
    calculator1Display.OSPRayScaleArray = 'Result'
    calculator1Display.OSPRayScaleFunction = 'PiecewiseFunction'
    calculator1Display.SelectOrientationVectors = 'U'
    calculator1Display.ScaleFactor = 0.10000629425048829
    calculator1Display.SelectScaleArray = 'Result'
    calculator1Display.GlyphType = 'Arrow'
    calculator1Display.GlyphTableIndexArray = 'Result'
    calculator1Display.GaussianRadius = 0.005000314712524414
    calculator1Display.SetScaleArray = ['POINTS', 'Result']
    calculator1Display.ScaleTransferFunction = 'PiecewiseFunction'
    calculator1Display.OpacityArray = ['POINTS', 'Result']
    calculator1Display.OpacityTransferFunction = 'PiecewiseFunction'
    calculator1Display.DataAxesGrid = 'GridAxesRepresentation'
    calculator1Display.PolarAxes = 'PolarAxesRepresentation'

    # init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
    calculator1Display.ScaleTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 14.734161467233076, 1.0, 0.5, 0.0]

    # init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
    calculator1Display.OpacityTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 14.734161467233076, 1.0, 0.5, 0.0]

    # hide data in view
    Hide(a0025foam, renderView1)

    # show color bar/color legend
    calculator1Display.SetScalarBarVisibility(renderView1, True)

    # update the view to ensure updated data information
    renderView1.Update()
    animationScene1 = GetAnimationScene()

    animationScene1.GoToLast()

    # get opacity transfer function/opacity map for 'Result'
    resultPWF = GetOpacityTransferFunction('Result')

    # create a new 'Integrate Variables'
    integrateVariables1 = IntegrateVariables(registrationName='IntegrateVariables1', Input=calculator1)

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

    # Properties modified on spreadSheetView1
    spreadSheetView1.ColumnToSort = 'Result'
    spreadSheetView1.InvertOrder = 1

    # Properties modified on spreadSheetView1
    spreadSheetView1.InvertOrder = 0

    # Properties modified on spreadSheetView1
    spreadSheetView1.InvertOrder = 1

    # export view
    ExportView(PATH+'FL/forces'+str(round(c,4)).replace('.','_')+'.csv', view=spreadSheetView1)

    # destroy spreadSheetView1
    Delete(spreadSheetView1)
    del spreadSheetView1

    # close an empty frame
    layout1.Collapse(2)

    # set active view
    SetActiveView(renderView1)

    # set active source
    SetActiveSource(calculator1)

    # destroy integrateVariables1
    Delete(integrateVariables1)
    del integrateVariables1

    # set active source
    SetActiveSource(a0025foam)

    # hide data in view
    Hide(calculator1, renderView1)

    # show data in view
    a0025foamDisplay = Show(a0025foam, renderView1, 'GeometryRepresentation')

    # show color bar/color legend
    a0025foamDisplay.SetScalarBarVisibility(renderView1, True)

    # destroy calculator1
    Delete(calculator1)
    del calculator1

for c in list(np.linspace(-L/2+dH/2,L/2-dH/2,20)):#list(np.linspace(-0.5,0.5,11)):
    # create a new 'Calculator'
    calculator1 = Calculator(registrationName='Calculator1', Input=a0025foam)
    calculator1.Function = ''
    # calculator1.ResultArrayName = 'Resul'+str(c)

    # Properties modified on calculator1
    C = c
    calculator1.Function = '(forces:moment_Z)*('+Cond+'<'+str(C+dH/2.)+')*('+Cond+'>'+str(C-dH/2.)+')'

    # show data in view
    calculator1Display = Show(calculator1, renderView1, 'GeometryRepresentation')

    # get color transfer function/color map for 'Result'
    resultLUT = GetColorTransferFunction('Result')

    # trace defaults for the display properties.
    calculator1Display.Representation = 'Surface'
    calculator1Display.ColorArrayName = ['POINTS', 'Result']
    calculator1Display.LookupTable = resultLUT
    calculator1Display.SelectTCoordArray = 'None'
    calculator1Display.SelectNormalArray = 'None'
    calculator1Display.SelectTangentArray = 'None'
    calculator1Display.OSPRayScaleArray = 'Result'
    calculator1Display.OSPRayScaleFunction = 'PiecewiseFunction'
    calculator1Display.SelectOrientationVectors = 'U'
    calculator1Display.ScaleFactor = 0.10000629425048829
    calculator1Display.SelectScaleArray = 'Result'
    calculator1Display.GlyphType = 'Arrow'
    calculator1Display.GlyphTableIndexArray = 'Result'
    calculator1Display.GaussianRadius = 0.005000314712524414
    calculator1Display.SetScaleArray = ['POINTS', 'Result']
    calculator1Display.ScaleTransferFunction = 'PiecewiseFunction'
    calculator1Display.OpacityArray = ['POINTS', 'Result']
    calculator1Display.OpacityTransferFunction = 'PiecewiseFunction'
    calculator1Display.DataAxesGrid = 'GridAxesRepresentation'
    calculator1Display.PolarAxes = 'PolarAxesRepresentation'

    # init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
    calculator1Display.ScaleTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 14.734161467233076, 1.0, 0.5, 0.0]

    # init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
    calculator1Display.OpacityTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 14.734161467233076, 1.0, 0.5, 0.0]

    # hide data in view
    Hide(a0025foam, renderView1)

    # show color bar/color legend
    calculator1Display.SetScalarBarVisibility(renderView1, True)

    # update the view to ensure updated data information
    renderView1.Update()
    animationScene1 = GetAnimationScene()

    animationScene1.GoToLast()

    # get opacity transfer function/opacity map for 'Result'
    resultPWF = GetOpacityTransferFunction('Result')

    # create a new 'Integrate Variables'
    integrateVariables1 = IntegrateVariables(registrationName='IntegrateVariables1', Input=calculator1)

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

    # Properties modified on spreadSheetView1
    spreadSheetView1.ColumnToSort = 'Result'
    spreadSheetView1.InvertOrder = 1

    # Properties modified on spreadSheetView1
    spreadSheetView1.InvertOrder = 0

    # Properties modified on spreadSheetView1
    spreadSheetView1.InvertOrder = 1

    # export view
    ExportView(PATH+'M/forces'+str(round(c,4)).replace('.','_')+'.csv', view=spreadSheetView1)

    # destroy spreadSheetView1
    Delete(spreadSheetView1)
    del spreadSheetView1

    # close an empty frame
    layout1.Collapse(2)

    # set active view
    SetActiveView(renderView1)

    # set active source
    SetActiveSource(calculator1)

    # destroy integrateVariables1
    Delete(integrateVariables1)
    del integrateVariables1

    # set active source
    SetActiveSource(a0025foam)

    # hide data in view
    Hide(calculator1, renderView1)

    # show data in view
    a0025foamDisplay = Show(a0025foam, renderView1, 'GeometryRepresentation')

    # show color bar/color legend
    a0025foamDisplay.SetScalarBarVisibility(renderView1, True)

    # destroy calculator1
    Delete(calculator1)
    del calculator1
#================================================================
# addendum: following script captures some of the application
# state to faithfully reproduce the visualization during playback
#================================================================

#--------------------------------
# saving layout sizes for layouts

# layout/tab size in pixels
layout1.SetSize(1153, 781)

#-----------------------------------
# saving camera placements for views

# current camera placement for renderView1
renderView1.CameraPosition = [12.0, 12.0, 14.429385609315649]
renderView1.CameraFocalPoint = [12.0, 12.0, 12.0]
renderView1.CameraParallelScale = 0.5196456723875056

#--------------------------------------------
# uncomment the following to render all views
# RenderAllViews()
# alternatively, if you want to write images, you can use SaveScreenshot(...).

paraview.simple._DisableFirstRenderCameraReset()

# get active source.
a0025foam = GetActiveSource()

# Properties modified on a0025foam
a0025foam.MeshRegions = ['air_to_spherePlaneFacesbot']
a0025foam.CellArrays = ['U', 'forces:force', 'forces:moment', 'p']

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
a0025foamDisplay.ScaleFactor = 0.10000629425048829
a0025foamDisplay.SelectScaleArray = 'p'
a0025foamDisplay.GlyphType = 'Arrow'
a0025foamDisplay.GlyphTableIndexArray = 'p'
a0025foamDisplay.GaussianRadius = 0.005000314712524414
a0025foamDisplay.SetScaleArray = ['POINTS', 'p']
a0025foamDisplay.ScaleTransferFunction = 'PiecewiseFunction'
a0025foamDisplay.OpacityArray = ['POINTS', 'p']
a0025foamDisplay.OpacityTransferFunction = 'PiecewiseFunction'
a0025foamDisplay.DataAxesGrid = 'GridAxesRepresentation'
a0025foamDisplay.PolarAxes = 'PolarAxesRepresentation'

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
a0025foamDisplay.ScaleTransferFunction.Points = [-105.09708404541016, 0.0, 0.5, 0.0, 103.13551330566406, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
a0025foamDisplay.OpacityTransferFunction.Points = [-105.09708404541016, 0.0, 0.5, 0.0, 103.13551330566406, 1.0, 0.5, 0.0]

# reset view to fit data
renderView1.ResetCamera()

# show color bar/color legend
a0025foamDisplay.SetScalarBarVisibility(renderView1, True)

# update the view to ensure updated data information
renderView1.Update()
animationScene1 = GetAnimationScene()

animationScene1.GoToLast()

# get opacity transfer function/opacity map for 'p'
pPWF = GetOpacityTransferFunction('p')

# set scalar coloring
ColorBy(a0025foamDisplay, ('POINTS', 'forces:force', 'Magnitude'))

# Hide the scalar bar for this color map if no visible data is colored by it.
HideScalarBarIfNotNeeded(pLUT, renderView1)

# rescale color and/or opacity maps used to include current data range
a0025foamDisplay.RescaleTransferFunctionToDataRange(True, False)

# show color bar/color legend
a0025foamDisplay.SetScalarBarVisibility(renderView1, True)

# get color transfer function/color map for 'forcesforce'
forcesforceLUT = GetColorTransferFunction('forcesforce')

# get opacity transfer function/opacity map for 'forcesforce'
forcesforcePWF = GetOpacityTransferFunction('forcesforce')







for c in [-L/2]:#list(np.linspace(-0.5,0.5,11)):
    # create a new 'Calculator'
    calculator1 = Calculator(registrationName='Calculator1', Input=a0025foam)
    calculator1.Function = ''
    # calculator1.ResultArrayName = 'Resul'+str(c)

    # Properties modified on calculator1
    C = c
    calculator1.Function = '(forces:force_X*'+str(ct)+'+'+str(st)+'*forces:force_Y)*('+Cond+'<'+str(C+dH/2.)+')*('+Cond+'>'+str(C-dH/2.)+')'

    # show data in view
    calculator1Display = Show(calculator1, renderView1, 'GeometryRepresentation')

    # get color transfer function/color map for 'Result'
    resultLUT = GetColorTransferFunction('Result')

    # trace defaults for the display properties.
    calculator1Display.Representation = 'Surface'
    calculator1Display.ColorArrayName = ['POINTS', 'Result']
    calculator1Display.LookupTable = resultLUT
    calculator1Display.SelectTCoordArray = 'None'
    calculator1Display.SelectNormalArray = 'None'
    calculator1Display.SelectTangentArray = 'None'
    calculator1Display.OSPRayScaleArray = 'Result'
    calculator1Display.OSPRayScaleFunction = 'PiecewiseFunction'
    calculator1Display.SelectOrientationVectors = 'U'
    calculator1Display.ScaleFactor = 0.10000629425048829
    calculator1Display.SelectScaleArray = 'Result'
    calculator1Display.GlyphType = 'Arrow'
    calculator1Display.GlyphTableIndexArray = 'Result'
    calculator1Display.GaussianRadius = 0.005000314712524414
    calculator1Display.SetScaleArray = ['POINTS', 'Result']
    calculator1Display.ScaleTransferFunction = 'PiecewiseFunction'
    calculator1Display.OpacityArray = ['POINTS', 'Result']
    calculator1Display.OpacityTransferFunction = 'PiecewiseFunction'
    calculator1Display.DataAxesGrid = 'GridAxesRepresentation'
    calculator1Display.PolarAxes = 'PolarAxesRepresentation'

    # init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
    calculator1Display.ScaleTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 14.734161467233076, 1.0, 0.5, 0.0]

    # init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
    calculator1Display.OpacityTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 14.734161467233076, 1.0, 0.5, 0.0]

    # hide data in view
    Hide(a0025foam, renderView1)

    # show color bar/color legend
    calculator1Display.SetScalarBarVisibility(renderView1, True)

    # update the view to ensure updated data information
    renderView1.Update()
    animationScene1 = GetAnimationScene()

    animationScene1.GoToLast()

    # get opacity transfer function/opacity map for 'Result'
    resultPWF = GetOpacityTransferFunction('Result')

    # create a new 'Integrate Variables'
    integrateVariables1 = IntegrateVariables(registrationName='IntegrateVariables1', Input=calculator1)

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

    # Properties modified on spreadSheetView1
    spreadSheetView1.ColumnToSort = 'Result'
    spreadSheetView1.InvertOrder = 1

    # Properties modified on spreadSheetView1
    spreadSheetView1.InvertOrder = 0

    # Properties modified on spreadSheetView1
    spreadSheetView1.InvertOrder = 1

    # export view
    ExportView(PATH+'FD/forces'+str(round(c,4)).replace('.','_')+'.csv', view=spreadSheetView1)


    # destroy spreadSheetView1
    Delete(spreadSheetView1)
    del spreadSheetView1

    # close an empty frame
    layout1.Collapse(2)

    # set active view
    SetActiveView(renderView1)

    # set active source
    SetActiveSource(calculator1)

    # destroy integrateVariables1
    Delete(integrateVariables1)
    del integrateVariables1

    # set active source
    SetActiveSource(a0025foam)

    # hide data in view
    Hide(calculator1, renderView1)

    # show data in view
    a0025foamDisplay = Show(a0025foam, renderView1, 'GeometryRepresentation')

    # show color bar/color legend
    a0025foamDisplay.SetScalarBarVisibility(renderView1, True)

    # destroy calculator1
    Delete(calculator1)
    del calculator1



for c in [-L/2]:#list(np.linspace(-0.5,0.5,11)):
    # create a new 'Calculator'
    calculator1 = Calculator(registrationName='Calculator1', Input=a0025foam)
    calculator1.Function = ''
    # calculator1.ResultArrayName = 'Resul'+str(c)

    # Properties modified on calculator1
    C = c
    calculator1.Function = '(-forces:force_X*'+str(st)+'+'+str(ct)+'*forces:force_Y)*('+Cond+'<'+str(C+dH/2.)+')*('+Cond+'>'+str(C-dH/2.)+')'

    # show data in view
    calculator1Display = Show(calculator1, renderView1, 'GeometryRepresentation')

    # get color transfer function/color map for 'Result'
    resultLUT = GetColorTransferFunction('Result')

    # trace defaults for the display properties.
    calculator1Display.Representation = 'Surface'
    calculator1Display.ColorArrayName = ['POINTS', 'Result']
    calculator1Display.LookupTable = resultLUT
    calculator1Display.SelectTCoordArray = 'None'
    calculator1Display.SelectNormalArray = 'None'
    calculator1Display.SelectTangentArray = 'None'
    calculator1Display.OSPRayScaleArray = 'Result'
    calculator1Display.OSPRayScaleFunction = 'PiecewiseFunction'
    calculator1Display.SelectOrientationVectors = 'U'
    calculator1Display.ScaleFactor = 0.10000629425048829
    calculator1Display.SelectScaleArray = 'Result'
    calculator1Display.GlyphType = 'Arrow'
    calculator1Display.GlyphTableIndexArray = 'Result'
    calculator1Display.GaussianRadius = 0.005000314712524414
    calculator1Display.SetScaleArray = ['POINTS', 'Result']
    calculator1Display.ScaleTransferFunction = 'PiecewiseFunction'
    calculator1Display.OpacityArray = ['POINTS', 'Result']
    calculator1Display.OpacityTransferFunction = 'PiecewiseFunction'
    calculator1Display.DataAxesGrid = 'GridAxesRepresentation'
    calculator1Display.PolarAxes = 'PolarAxesRepresentation'

    # init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
    calculator1Display.ScaleTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 14.734161467233076, 1.0, 0.5, 0.0]

    # init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
    calculator1Display.OpacityTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 14.734161467233076, 1.0, 0.5, 0.0]

    # hide data in view
    Hide(a0025foam, renderView1)

    # show color bar/color legend
    calculator1Display.SetScalarBarVisibility(renderView1, True)

    # update the view to ensure updated data information
    renderView1.Update()
    animationScene1 = GetAnimationScene()

    animationScene1.GoToLast()

    # get opacity transfer function/opacity map for 'Result'
    resultPWF = GetOpacityTransferFunction('Result')

    # create a new 'Integrate Variables'
    integrateVariables1 = IntegrateVariables(registrationName='IntegrateVariables1', Input=calculator1)

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

    # Properties modified on spreadSheetView1
    spreadSheetView1.ColumnToSort = 'Result'
    spreadSheetView1.InvertOrder = 1

    # Properties modified on spreadSheetView1
    spreadSheetView1.InvertOrder = 0

    # Properties modified on spreadSheetView1
    spreadSheetView1.InvertOrder = 1

    # export view
    ExportView(PATH+'FL/forces'+str(round(c,4)).replace('.','_')+'.csv', view=spreadSheetView1)

    # destroy spreadSheetView1
    Delete(spreadSheetView1)
    del spreadSheetView1

    # close an empty frame
    layout1.Collapse(2)

    # set active view
    SetActiveView(renderView1)

    # set active source
    SetActiveSource(calculator1)

    # destroy integrateVariables1
    Delete(integrateVariables1)
    del integrateVariables1

    # set active source
    SetActiveSource(a0025foam)

    # hide data in view
    Hide(calculator1, renderView1)

    # show data in view
    a0025foamDisplay = Show(a0025foam, renderView1, 'GeometryRepresentation')

    # show color bar/color legend
    a0025foamDisplay.SetScalarBarVisibility(renderView1, True)

    # destroy calculator1
    Delete(calculator1)
    del calculator1

for c in [-L/2]:#list(np.linspace(-0.5,0.5,11)):
    # create a new 'Calculator'
    calculator1 = Calculator(registrationName='Calculator1', Input=a0025foam)
    calculator1.Function = ''
    # calculator1.ResultArrayName = 'Resul'+str(c)

    # Properties modified on calculator1
    C = c
    calculator1.Function = '(forces:moment_Z)*('+Cond+'<'+str(C+dH/2.)+')*('+Cond+'>'+str(C-dH/2.)+')'

    # show data in view
    calculator1Display = Show(calculator1, renderView1, 'GeometryRepresentation')

    # get color transfer function/color map for 'Result'
    resultLUT = GetColorTransferFunction('Result')

    # trace defaults for the display properties.
    calculator1Display.Representation = 'Surface'
    calculator1Display.ColorArrayName = ['POINTS', 'Result']
    calculator1Display.LookupTable = resultLUT
    calculator1Display.SelectTCoordArray = 'None'
    calculator1Display.SelectNormalArray = 'None'
    calculator1Display.SelectTangentArray = 'None'
    calculator1Display.OSPRayScaleArray = 'Result'
    calculator1Display.OSPRayScaleFunction = 'PiecewiseFunction'
    calculator1Display.SelectOrientationVectors = 'U'
    calculator1Display.ScaleFactor = 0.10000629425048829
    calculator1Display.SelectScaleArray = 'Result'
    calculator1Display.GlyphType = 'Arrow'
    calculator1Display.GlyphTableIndexArray = 'Result'
    calculator1Display.GaussianRadius = 0.005000314712524414
    calculator1Display.SetScaleArray = ['POINTS', 'Result']
    calculator1Display.ScaleTransferFunction = 'PiecewiseFunction'
    calculator1Display.OpacityArray = ['POINTS', 'Result']
    calculator1Display.OpacityTransferFunction = 'PiecewiseFunction'
    calculator1Display.DataAxesGrid = 'GridAxesRepresentation'
    calculator1Display.PolarAxes = 'PolarAxesRepresentation'

    # init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
    calculator1Display.ScaleTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 14.734161467233076, 1.0, 0.5, 0.0]

    # init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
    calculator1Display.OpacityTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 14.734161467233076, 1.0, 0.5, 0.0]

    # hide data in view
    Hide(a0025foam, renderView1)

    # show color bar/color legend
    calculator1Display.SetScalarBarVisibility(renderView1, True)

    # update the view to ensure updated data information
    renderView1.Update()
    animationScene1 = GetAnimationScene()

    animationScene1.GoToLast()

    # get opacity transfer function/opacity map for 'Result'
    resultPWF = GetOpacityTransferFunction('Result')

    # create a new 'Integrate Variables'
    integrateVariables1 = IntegrateVariables(registrationName='IntegrateVariables1', Input=calculator1)

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

    # Properties modified on spreadSheetView1
    spreadSheetView1.ColumnToSort = 'Result'
    spreadSheetView1.InvertOrder = 1

    # Properties modified on spreadSheetView1
    spreadSheetView1.InvertOrder = 0

    # Properties modified on spreadSheetView1
    spreadSheetView1.InvertOrder = 1

    # export view
    ExportView(PATH+'M/forces'+str(round(c,4)).replace('.','_')+'.csv', view=spreadSheetView1)

    # destroy spreadSheetView1
    Delete(spreadSheetView1)
    del spreadSheetView1

    # close an empty frame
    layout1.Collapse(2)

    # set active view
    SetActiveView(renderView1)

    # set active source
    SetActiveSource(calculator1)

    # destroy integrateVariables1
    Delete(integrateVariables1)
    del integrateVariables1

    # set active source
    SetActiveSource(a0025foam)

    # hide data in view
    Hide(calculator1, renderView1)

    # show data in view
    a0025foamDisplay = Show(a0025foam, renderView1, 'GeometryRepresentation')

    # show color bar/color legend
    a0025foamDisplay.SetScalarBarVisibility(renderView1, True)

    # destroy calculator1
    Delete(calculator1)
    del calculator1
#================================================================
# addendum: following script captures some of the application
# state to faithfully reproduce the visualization during playback
#================================================================

#--------------------------------
# saving layout sizes for layouts

# layout/tab size in pixels
layout1.SetSize(1153, 781)

#-----------------------------------
# saving camera placements for views

# current camera placement for renderView1
renderView1.CameraPosition = [12.0, 12.0, 14.429385609315649]
renderView1.CameraFocalPoint = [12.0, 12.0, 12.0]
renderView1.CameraParallelScale = 0.5196456723875056

#--------------------------------------------
# uncomment the following to render all views
# RenderAllViews()
# alternatively, if you want to write images, you can use SaveScreenshot(...).

paraview.simple._DisableFirstRenderCameraReset()

# get active source.
a0025foam = GetActiveSource()

# Properties modified on a0025foam
a0025foam.MeshRegions = ['air_to_spherePlaneFacestop']
a0025foam.CellArrays = ['U', 'forces:force', 'forces:moment', 'p']

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
a0025foamDisplay.ScaleFactor = 0.10000629425048829
a0025foamDisplay.SelectScaleArray = 'p'
a0025foamDisplay.GlyphType = 'Arrow'
a0025foamDisplay.GlyphTableIndexArray = 'p'
a0025foamDisplay.GaussianRadius = 0.005000314712524414
a0025foamDisplay.SetScaleArray = ['POINTS', 'p']
a0025foamDisplay.ScaleTransferFunction = 'PiecewiseFunction'
a0025foamDisplay.OpacityArray = ['POINTS', 'p']
a0025foamDisplay.OpacityTransferFunction = 'PiecewiseFunction'
a0025foamDisplay.DataAxesGrid = 'GridAxesRepresentation'
a0025foamDisplay.PolarAxes = 'PolarAxesRepresentation'

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
a0025foamDisplay.ScaleTransferFunction.Points = [-105.09708404541016, 0.0, 0.5, 0.0, 103.13551330566406, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
a0025foamDisplay.OpacityTransferFunction.Points = [-105.09708404541016, 0.0, 0.5, 0.0, 103.13551330566406, 1.0, 0.5, 0.0]

# reset view to fit data
renderView1.ResetCamera()

# show color bar/color legend
a0025foamDisplay.SetScalarBarVisibility(renderView1, True)

# update the view to ensure updated data information
renderView1.Update()
animationScene1 = GetAnimationScene()

animationScene1.GoToLast()

# get opacity transfer function/opacity map for 'p'
pPWF = GetOpacityTransferFunction('p')

# set scalar coloring
ColorBy(a0025foamDisplay, ('POINTS', 'forces:force', 'Magnitude'))

# Hide the scalar bar for this color map if no visible data is colored by it.
HideScalarBarIfNotNeeded(pLUT, renderView1)

# rescale color and/or opacity maps used to include current data range
a0025foamDisplay.RescaleTransferFunctionToDataRange(True, False)

# show color bar/color legend
a0025foamDisplay.SetScalarBarVisibility(renderView1, True)

# get color transfer function/color map for 'forcesforce'
forcesforceLUT = GetColorTransferFunction('forcesforce')

# get opacity transfer function/opacity map for 'forcesforce'
forcesforcePWF = GetOpacityTransferFunction('forcesforce')

####################
####################
####################
####################
####################
####################





for c in [L/2]:#list(np.linspace(-0.5,0.5,11)):
    # create a new 'Calculator'
    calculator1 = Calculator(registrationName='Calculator1', Input=a0025foam)
    calculator1.Function = ''
    # calculator1.ResultArrayName = 'Resul'+str(c)

    # Properties modified on calculator1
    C = c
    calculator1.Function = '(forces:force_X*'+str(ct)+'+'+str(st)+'*forces:force_Y)*('+Cond+'<'+str(C+dH/2.)+')*('+Cond+'>'+str(C-dH/2.)+')'

    # show data in view
    calculator1Display = Show(calculator1, renderView1, 'GeometryRepresentation')

    # get color transfer function/color map for 'Result'
    resultLUT = GetColorTransferFunction('Result')

    # trace defaults for the display properties.
    calculator1Display.Representation = 'Surface'
    calculator1Display.ColorArrayName = ['POINTS', 'Result']
    calculator1Display.LookupTable = resultLUT
    calculator1Display.SelectTCoordArray = 'None'
    calculator1Display.SelectNormalArray = 'None'
    calculator1Display.SelectTangentArray = 'None'
    calculator1Display.OSPRayScaleArray = 'Result'
    calculator1Display.OSPRayScaleFunction = 'PiecewiseFunction'
    calculator1Display.SelectOrientationVectors = 'U'
    calculator1Display.ScaleFactor = 0.10000629425048829
    calculator1Display.SelectScaleArray = 'Result'
    calculator1Display.GlyphType = 'Arrow'
    calculator1Display.GlyphTableIndexArray = 'Result'
    calculator1Display.GaussianRadius = 0.005000314712524414
    calculator1Display.SetScaleArray = ['POINTS', 'Result']
    calculator1Display.ScaleTransferFunction = 'PiecewiseFunction'
    calculator1Display.OpacityArray = ['POINTS', 'Result']
    calculator1Display.OpacityTransferFunction = 'PiecewiseFunction'
    calculator1Display.DataAxesGrid = 'GridAxesRepresentation'
    calculator1Display.PolarAxes = 'PolarAxesRepresentation'

    # init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
    calculator1Display.ScaleTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 14.734161467233076, 1.0, 0.5, 0.0]

    # init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
    calculator1Display.OpacityTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 14.734161467233076, 1.0, 0.5, 0.0]

    # hide data in view
    Hide(a0025foam, renderView1)

    # show color bar/color legend
    calculator1Display.SetScalarBarVisibility(renderView1, True)

    # update the view to ensure updated data information
    renderView1.Update()
    animationScene1 = GetAnimationScene()

    animationScene1.GoToLast()

    # get opacity transfer function/opacity map for 'Result'
    resultPWF = GetOpacityTransferFunction('Result')

    # create a new 'Integrate Variables'
    integrateVariables1 = IntegrateVariables(registrationName='IntegrateVariables1', Input=calculator1)

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

    # Properties modified on spreadSheetView1
    spreadSheetView1.ColumnToSort = 'Result'
    spreadSheetView1.InvertOrder = 1

    # Properties modified on spreadSheetView1
    spreadSheetView1.InvertOrder = 0

    # Properties modified on spreadSheetView1
    spreadSheetView1.InvertOrder = 1

    # export view
    ExportView(PATH+'FD/forces'+str(round(c,4)).replace('.','_')+'.csv', view=spreadSheetView1)


    # destroy spreadSheetView1
    Delete(spreadSheetView1)
    del spreadSheetView1

    # close an empty frame
    layout1.Collapse(2)

    # set active view
    SetActiveView(renderView1)

    # set active source
    SetActiveSource(calculator1)

    # destroy integrateVariables1
    Delete(integrateVariables1)
    del integrateVariables1

    # set active source
    SetActiveSource(a0025foam)

    # hide data in view
    Hide(calculator1, renderView1)

    # show data in view
    a0025foamDisplay = Show(a0025foam, renderView1, 'GeometryRepresentation')

    # show color bar/color legend
    a0025foamDisplay.SetScalarBarVisibility(renderView1, True)

    # destroy calculator1
    Delete(calculator1)
    del calculator1




for c in [L/2]:#list(np.linspace(-0.5,0.5,11)):
    # create a new 'Calculator'
    calculator1 = Calculator(registrationName='Calculator1', Input=a0025foam)
    calculator1.Function = ''
    # calculator1.ResultArrayName = 'Resul'+str(c)

    # Properties modified on calculator1
    C = c
    calculator1.Function = '(-forces:force_X*'+str(st)+'+'+str(ct)+'*forces:force_Y)*('+Cond+'<'+str(C+dH/2.)+')*('+Cond+'>'+str(C-dH/2.)+')'

    # show data in view
    calculator1Display = Show(calculator1, renderView1, 'GeometryRepresentation')

    # get color transfer function/color map for 'Result'
    resultLUT = GetColorTransferFunction('Result')

    # trace defaults for the display properties.
    calculator1Display.Representation = 'Surface'
    calculator1Display.ColorArrayName = ['POINTS', 'Result']
    calculator1Display.LookupTable = resultLUT
    calculator1Display.SelectTCoordArray = 'None'
    calculator1Display.SelectNormalArray = 'None'
    calculator1Display.SelectTangentArray = 'None'
    calculator1Display.OSPRayScaleArray = 'Result'
    calculator1Display.OSPRayScaleFunction = 'PiecewiseFunction'
    calculator1Display.SelectOrientationVectors = 'U'
    calculator1Display.ScaleFactor = 0.10000629425048829
    calculator1Display.SelectScaleArray = 'Result'
    calculator1Display.GlyphType = 'Arrow'
    calculator1Display.GlyphTableIndexArray = 'Result'
    calculator1Display.GaussianRadius = 0.005000314712524414
    calculator1Display.SetScaleArray = ['POINTS', 'Result']
    calculator1Display.ScaleTransferFunction = 'PiecewiseFunction'
    calculator1Display.OpacityArray = ['POINTS', 'Result']
    calculator1Display.OpacityTransferFunction = 'PiecewiseFunction'
    calculator1Display.DataAxesGrid = 'GridAxesRepresentation'
    calculator1Display.PolarAxes = 'PolarAxesRepresentation'

    # init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
    calculator1Display.ScaleTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 14.734161467233076, 1.0, 0.5, 0.0]

    # init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
    calculator1Display.OpacityTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 14.734161467233076, 1.0, 0.5, 0.0]

    # hide data in view
    Hide(a0025foam, renderView1)

    # show color bar/color legend
    calculator1Display.SetScalarBarVisibility(renderView1, True)

    # update the view to ensure updated data information
    renderView1.Update()
    animationScene1 = GetAnimationScene()

    animationScene1.GoToLast()

    # get opacity transfer function/opacity map for 'Result'
    resultPWF = GetOpacityTransferFunction('Result')

    # create a new 'Integrate Variables'
    integrateVariables1 = IntegrateVariables(registrationName='IntegrateVariables1', Input=calculator1)

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

    # Properties modified on spreadSheetView1
    spreadSheetView1.ColumnToSort = 'Result'
    spreadSheetView1.InvertOrder = 1

    # Properties modified on spreadSheetView1
    spreadSheetView1.InvertOrder = 0

    # Properties modified on spreadSheetView1
    spreadSheetView1.InvertOrder = 1

    # export view
    ExportView(PATH+'FL/forces'+str(round(c,4)).replace('.','_')+'.csv', view=spreadSheetView1)
    print('j ai exporter le sbailss')
    # destroy spreadSheetView1
    Delete(spreadSheetView1)
    del spreadSheetView1

    # close an empty frame
    layout1.Collapse(2)

    # set active view
    SetActiveView(renderView1)

    # set active source
    SetActiveSource(calculator1)

    # destroy integrateVariables1
    Delete(integrateVariables1)
    del integrateVariables1

    # set active source
    SetActiveSource(a0025foam)

    # hide data in view
    Hide(calculator1, renderView1)

    # show data in view
    a0025foamDisplay = Show(a0025foam, renderView1, 'GeometryRepresentation')

    # show color bar/color legend
    a0025foamDisplay.SetScalarBarVisibility(renderView1, True)

    # destroy calculator1
    Delete(calculator1)
    del calculator1

for c in [L/2]:#list(np.linspace(-0.5,0.5,11)):
    # create a new 'Calculator'
    calculator1 = Calculator(registrationName='Calculator1', Input=a0025foam)
    calculator1.Function = ''
    # calculator1.ResultArrayName = 'Resul'+str(c)

    # Properties modified on calculator1
    C = c
    calculator1.Function = '(forces:moment_Z)*('+Cond+'<'+str(C+dH/2.)+')*('+Cond+'>'+str(C-dH/2.)+')'

    # show data in view
    calculator1Display = Show(calculator1, renderView1, 'GeometryRepresentation')

    # get color transfer function/color map for 'Result'
    resultLUT = GetColorTransferFunction('Result')

    # trace defaults for the display properties.
    calculator1Display.Representation = 'Surface'
    calculator1Display.ColorArrayName = ['POINTS', 'Result']
    calculator1Display.LookupTable = resultLUT
    calculator1Display.SelectTCoordArray = 'None'
    calculator1Display.SelectNormalArray = 'None'
    calculator1Display.SelectTangentArray = 'None'
    calculator1Display.OSPRayScaleArray = 'Result'
    calculator1Display.OSPRayScaleFunction = 'PiecewiseFunction'
    calculator1Display.SelectOrientationVectors = 'U'
    calculator1Display.ScaleFactor = 0.10000629425048829
    calculator1Display.SelectScaleArray = 'Result'
    calculator1Display.GlyphType = 'Arrow'
    calculator1Display.GlyphTableIndexArray = 'Result'
    calculator1Display.GaussianRadius = 0.005000314712524414
    calculator1Display.SetScaleArray = ['POINTS', 'Result']
    calculator1Display.ScaleTransferFunction = 'PiecewiseFunction'
    calculator1Display.OpacityArray = ['POINTS', 'Result']
    calculator1Display.OpacityTransferFunction = 'PiecewiseFunction'
    calculator1Display.DataAxesGrid = 'GridAxesRepresentation'
    calculator1Display.PolarAxes = 'PolarAxesRepresentation'

    # init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
    calculator1Display.ScaleTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 14.734161467233076, 1.0, 0.5, 0.0]

    # init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
    calculator1Display.OpacityTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 14.734161467233076, 1.0, 0.5, 0.0]

    # hide data in view
    Hide(a0025foam, renderView1)

    # show color bar/color legend
    calculator1Display.SetScalarBarVisibility(renderView1, True)

    # update the view to ensure updated data information
    renderView1.Update()
    animationScene1 = GetAnimationScene()

    animationScene1.GoToLast()

    # get opacity transfer function/opacity map for 'Result'
    resultPWF = GetOpacityTransferFunction('Result')

    # create a new 'Integrate Variables'
    integrateVariables1 = IntegrateVariables(registrationName='IntegrateVariables1', Input=calculator1)

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

    # Properties modified on spreadSheetView1
    spreadSheetView1.ColumnToSort = 'Result'
    spreadSheetView1.InvertOrder = 1

    # Properties modified on spreadSheetView1
    spreadSheetView1.InvertOrder = 0

    # Properties modified on spreadSheetView1
    spreadSheetView1.InvertOrder = 1

    # export view
    ExportView(PATH+'M/forces'+str(round(c,4)).replace('.','_')+'.csv', view=spreadSheetView1)
    print('j ai exporter le sbailss')
    # destroy spreadSheetView1
    Delete(spreadSheetView1)
    del spreadSheetView1

    # close an empty frame
    layout1.Collapse(2)

    # set active view
    SetActiveView(renderView1)

    # set active source
    SetActiveSource(calculator1)

    # destroy integrateVariables1
    Delete(integrateVariables1)
    del integrateVariables1

    # set active source
    SetActiveSource(a0025foam)

    # hide data in view
    Hide(calculator1, renderView1)

    # show data in view
    a0025foamDisplay = Show(a0025foam, renderView1, 'GeometryRepresentation')

    # show color bar/color legend
    a0025foamDisplay.SetScalarBarVisibility(renderView1, True)

    # destroy calculator1
    Delete(calculator1)
    del calculator1
#================================================================
# addendum: following script captures some of the application
# state to faithfully reproduce the visualization during playback
#================================================================

#--------------------------------
# saving layout sizes for layouts

# layout/tab size in pixels
layout1.SetSize(1153, 781)

#-----------------------------------
# saving camera placements for views

# current camera placement for renderView1
renderView1.CameraPosition = [12.0, 12.0, 14.429385609315649]
renderView1.CameraFocalPoint = [12.0, 12.0, 12.0]
renderView1.CameraParallelScale = 0.5196456723875056

#--------------------------------------------
# uncomment the following to render all views
# RenderAllViews()
# alternatively, if you want to write images, you can use SaveScreenshot(...).