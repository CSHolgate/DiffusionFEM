#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# create a new 'PVD Reader'
solutionpvd = PVDReader(FileName='/Users/Collin/Fenics/Diffusion/ExperimentalRecreation/1D/data1D_ER/solution.pvd')

# get animation scene
animationScene1 = GetAnimationScene()

# update animation scene based on data timesteps
animationScene1.UpdateAnimationUsingDataTimeSteps()

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')
# uncomment following to set a specific view size
# renderView1.ViewSize = [1011, 550]

# get color transfer function/color map for 'f_10'
f_10LUT = GetColorTransferFunction('f_10')

# show data in view
solutionpvdDisplay = Show(solutionpvd, renderView1)
# trace defaults for the display properties.
solutionpvdDisplay.Representation = 'Surface'
solutionpvdDisplay.ColorArrayName = ['POINTS', 'f_10']
solutionpvdDisplay.LookupTable = f_10LUT
solutionpvdDisplay.OSPRayScaleArray = 'f_10'
solutionpvdDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
solutionpvdDisplay.SelectOrientationVectors = 'None'
solutionpvdDisplay.ScaleFactor = 0.2
solutionpvdDisplay.SelectScaleArray = 'f_10'
solutionpvdDisplay.GlyphType = 'Arrow'
solutionpvdDisplay.PolarAxes = 'PolarAxesRepresentation'
solutionpvdDisplay.ScalarOpacityUnitDistance = 0.20000000000000004
solutionpvdDisplay.GaussianRadius = 0.1
solutionpvdDisplay.SetScaleArray = ['POINTS', 'f_10']
solutionpvdDisplay.ScaleTransferFunction = 'PiecewiseFunction'
solutionpvdDisplay.OpacityArray = ['POINTS', 'f_10']
solutionpvdDisplay.OpacityTransferFunction = 'PiecewiseFunction'

# reset view to fit data
renderView1.ResetCamera()

# show color bar/color legend
solutionpvdDisplay.SetScalarBarVisibility(renderView1, True)

animationScene1.GoToLast()

# save data
SaveData('/Users/Collin/Fenics/Diffusion/ExperimentalRecreation/1D/data1D_ER/Fenics.csv', proxy=solutionpvd)

#### saving camera placements for all active views

# current camera placement for renderView1
renderView1.InteractionMode = '2D'
renderView1.CameraPosition = [1.0, 3.8637033051562737, 0.0]
renderView1.CameraFocalPoint = [1.0, 0.0, 0.0]
renderView1.CameraViewUp = [1.0, 0.0, 0.0]
renderView1.CameraParallelScale = 1.0

#### uncomment the following to render all views
# RenderAllViews()
# alternatively, if you want to write images, you can use SaveScreenshot(...).