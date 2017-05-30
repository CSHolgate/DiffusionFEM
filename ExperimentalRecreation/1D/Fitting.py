from __future__ import print_function
import numpy
import scipy
import fenics
from scipy.optimize import curve_fit
from fenics import *
from scipy.interpolate import interp1d

try:
	nodes = int(input('How many nodes to run? (Default value is 100) '))
except SyntaxError:
	nodes = 100

def func(x, Co, D):
	T = 3480.0
	num_steps = 300
	dt = T/num_steps

	diff = D

	# Print the present guesses for Co and D
	print(diff)
	print(Co)

	# Create the mesh and define function space
	mesh = IntervalMesh(nodes, 0, 2)
	V = FunctionSpace(mesh, 'P',1)

	# Define the tolerance of the boundary condition
	tol = 1.0E-14

	# Define the boundary condition
	u_D = Expression('Coo * (1- exp(-0.49 * (t/60)))',degree=2, t=0, Coo=Co)

	def boundary_D(x, on_boundary):
		if on_boundary:
			if near(x[0], 0, tol):
				return True
			else:
				return False
		else:
			return False

	bc = DirichletBC(V, u_D, boundary_D)

	# Define the initial value
	u_n = interpolate(u_D, V)

	# Define the variational problem
	u = TrialFunction(V)
	v = TestFunction(V)
	f = Constant(0.0)

	a = (u*v + dt*diff*dot(grad(u), grad(v)))*dx
	L = (u_n + dt*f)*v*dx

	# Create VTK file for saving solution
	#vtkfile = File('dataDiffFConc/solution.pvd')

	# Time stepping
	u = Function(V)
	t = 0.0

	for n in range(num_steps):
		t += dt
		u_D.t = t

		# Compute Solution
		solve(a == L, u, bc)

		# Save to file and plot solution
		#vtkfile << (u,t)

		# Compute Error at vertices
		# Not done here

		# Print time
		print('t = %.2f' % t)

		# Update previous solution
		u_n.assign(u)

	print('\n' * 10)

	# Run the fitting methods
	interfun = interp1d(numpy.linspace(2,0,num = nodes+1), u.vector().array())
	return interfun(x)

# Data PX-YSZ, Kramer CMAS, 1300C, 60 Minutes data for Zr
xdata = [18.4631,24.1199,29.7768,36.1799,41.8368,48.2399,53.8967,60.2999,65.9567,72.3598,78.017,84.42,90.077,95.734,102.137,108.54,114.197,119.854,126.257,132.66,138.317,143.973,151.045,156.701,162.358,168.015,174.418,180.821,186.478,192.135,198.538,204.941,210.598,216.255,222.658,228.315,234.718,240.375,246.778,252.435,258.838,264.495,270.898,276.555,282.212,289.283,294.94,300.597,306.254,313.325,318.981,324.638,330.295,336.698,343.101,348.758,354.415,360.818,367.221,372.878,378.535,384.938,390.595,396.998,402.655,409.058,414.715,421.118,426.775,433.178,438.835,445.238,451.641,457.298,462.955,468.612,475.683,481.34,486.997,492.654]
xdata[:] = [x / 1000 for x in xdata]
ydata = [2.59114926,2.48375488,2.432533019,2.240920246,2.232590513,2.161828643,2.059430461,1.996965007,1.887728518,1.907844208,1.808829964,1.57142775,1.701358373,1.605451805,1.456515395,1.461104922,1.352972172,1.302182066,1.243359915,1.066671415,1.159614932,1.126754173,1.070932217,1.032455722,1.045294806,0.96000336,0.893994653,0.903897486,0.830675659,0.774701278,0.718758126,0.716032941,0.713599868,0.630989886,0.648870892,0.567879499,0.541446439,0.548667775,0.468353856,0.495793694,0.450233883,0.435563393,0.403205311,0.427295594,0.359591388,-0.582969908,0.335402222,0.309952966,0.300129552,0.293129426,0.293583798,0.239513564,0.255303141,0.220204387,0.225505408,0.227246057,0.154325388,0.199752451,0.169280343,0.176321289,0.167477769,0.145899825,0.125536605,0.12135882,0.114396878,0.11778856,0.118246595,0.102681041,0.085488536,0.070239495,0.093348714,0.088575438,0.077767637,0.078043344,0.071393113,0.041080937,0.058835859,0.048843231,0.051745413,0.06697217]

best_fit_data,fit_errors = curve_fit(func, xdata, ydata,[2.8,4E-6])

print(best_fit_data)
print(fit_errors)
print()
perr = numpy.sqrt(numpy.diag(fit_errors))
print(perr)

a1,a2 = best_fit_data
print("The best fit parameter values are: ")
print("Co = {}".format(a1))
print("Do = {}".format(a2))









