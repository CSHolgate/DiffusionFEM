#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# create a new 'PVD Reader'
solutionpvd = PVDReader(FileName='/Users/Collin/Fenics/Diffusion/NonlinearD/dataDiffFConc/solution.pvd')

# get animation scene
animationScene1 = GetAnimationScene()

# update animation scene based on data timesteps
animationScene1.UpdateAnimationUsingDataTimeSteps()

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')
# uncomment following to set a specific view size
# renderView1.ViewSize = [1011, 550]

# get color transfer function/color map for 'f_9'
f_9LUT = GetColorTransferFunction('f_9')

# show data in view
solutionpvdDisplay = Show(solutionpvd, renderView1)
# trace defaults for the display properties.
solutionpvdDisplay.Representation = 'Surface'
solutionpvdDisplay.ColorArrayName = ['POINTS', 'f_9']
solutionpvdDisplay.LookupTable = f_9LUT
solutionpvdDisplay.OSPRayScaleArray = 'f_9'
solutionpvdDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
solutionpvdDisplay.SelectOrientationVectors = 'None'
solutionpvdDisplay.ScaleFactor = 0.2
solutionpvdDisplay.SelectScaleArray = 'f_9'
solutionpvdDisplay.GlyphType = 'Arrow'
solutionpvdDisplay.PolarAxes = 'PolarAxesRepresentation'
solutionpvdDisplay.ScalarOpacityUnitDistance = 0.10707485778544379
solutionpvdDisplay.GaussianRadius = 0.1
solutionpvdDisplay.SetScaleArray = ['POINTS', 'f_9']
solutionpvdDisplay.ScaleTransferFunction = 'PiecewiseFunction'
solutionpvdDisplay.OpacityArray = ['POINTS', 'f_9']
solutionpvdDisplay.OpacityTransferFunction = 'PiecewiseFunction'

# reset view to fit data
renderView1.ResetCamera()

#changing interaction mode based on data extents
renderView1.CameraPosition = [0.0, 1.0, 10000.0]

# show color bar/color legend
solutionpvdDisplay.SetScalarBarVisibility(renderView1, True)

#animationScene1.Play()

animationScene1.GoToLast()

animationScene1.GoToPrevious()

animationScene1.GoToPrevious()

animationScene1.GoToPrevious()

animationScene1.GoToPrevious()

animationScene1.GoToPrevious()

animationScene1.GoToPrevious()

animationScene1.GoToPrevious()

animationScene1.GoToPrevious()

animationScene1.GoToPrevious()

animationScene1.GoToPrevious()

# save data
SaveData('/Users/Collin/Fenics/Diffusion/NonlinearD/dataDiffFConc/fenics.csv', proxy=solutionpvd)

#### saving camera placements for all active views

# current camera placement for renderView1
renderView1.InteractionMode = '2D'
renderView1.CameraPosition = [0.0, 1.0, 10000.0]
renderView1.CameraFocalPoint = [0.0, 1.0, 0.0]
renderView1.CameraParallelScale = 1.4142135623730951

#### uncomment the following to render all views
# RenderAllViews()
# alternatively, if you want to write images, you can use SaveScreenshot(...).