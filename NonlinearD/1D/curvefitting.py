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


def func(x, a, Co, do):
	# Define Time Conditions
	T = 3480.0
	num_steps = 300
	dt = T/num_steps

	# Define the non-linear coefficient Diff
	def Diff(u):
		"Return nonlinear coefficient"
		#return 5.0E-6 * (1.0 - 0.36*u + ((0.36*u)**2.0)/2.0)
		return do * fenics.exp(a * u)

	# Display the present parameter guesses of a, Co, and do
	print(a)
	print(Co)
	print(do)


	# Create mesh and define function space
	mesh = IntervalMesh(nodes, 0, 2)
	V = FunctionSpace(mesh, 'P', 1)

	# Define the tolerance of the boundary condition
	tol = 1.0E-14

	# Define the boundary condition
	# Including mixed Neumann (0 at boundaries) and Dirichlet conditions
	u_D = Expression('Coo * (1 - exp(-30000 * t))', degree=2, t=0, Coo=Co) # DEGREE?

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
	u = Function(V) # Note: NOT TRIALFUNCTION!
	v = TestFunction(V)
	f = Constant(0.0)

	#a = (u*v + dt*Diff(u)*dot(grad(u), grad(v)))*dx
	#L = (u_n + dt*f)*v*dx
	F = (u*v + dt*Diff(u)*dot(grad(u),grad(v)))*dx - (u_n + dt*f)*v*dx

	# Create VTK file for saving solution
	#vtkfile = File('dataDiffFConc/solution.pvd')

	# Time stepping
	#u = Function(V)
	t = 0.0

	for n in range(num_steps):

		# Update current time
		t += dt
		u_D.t = t

		# Print time
		print('t= %2.f' % t)

		# Compute solution
		solve(F == 0, u, bc)

		# Save to file and plot solution
		#vtkfile << (u,t)

		# Update previous solution
		u_n.assign(u)
	
	print("""










		""")
	#print(u.vector().array().tolist())
	#print(numpy.linspace(2,0,num=101).tolist())
	interfun = interp1d(numpy.linspace(2,0,num=nodes+1),u.vector().array())
	return interfun(x)

	
# PX-YSZ, KRAMER CMAS, 1300C, 60 MINUTES DATA FOR ZR

#xdata = [18.4631,24.1199,29.7768,36.1799,41.8368,48.2399,53.8967,60.2999,65.9567,72.3598,78.017,84.42,90.077,95.734,102.137,108.54,114.197,119.854,126.257,132.66,138.317,143.973,151.045,156.701,162.358,168.015,174.418,180.821,186.478,192.135,198.538,204.941,210.598,216.255,222.658,228.315,234.718,240.375,246.778,252.435,258.838,264.495,270.898,276.555,282.212,289.283,294.94,300.597,306.254,313.325,318.981,324.638,330.295,336.698,343.101,348.758,354.415,360.818,367.221,372.878,378.535,384.938,390.595,396.998,402.655,409.058,414.715,421.118,426.775,433.178,438.835,445.238,451.641,457.298,462.955,468.612,475.683,481.34,486.997,492.654]
#xdata[:] = [x / 1000 for x in xdata]
#ydata = [2.59114926,2.48375488,2.432533019,2.240920246,2.232590513,2.161828643,2.059430461,1.996965007,1.887728518,1.907844208,1.808829964,1.57142775,1.701358373,1.605451805,1.456515395,1.461104922,1.352972172,1.302182066,1.243359915,1.066671415,1.159614932,1.126754173,1.070932217,1.032455722,1.045294806,0.96000336,0.893994653,0.903897486,0.830675659,0.774701278,0.718758126,0.716032941,0.713599868,0.630989886,0.648870892,0.567879499,0.541446439,0.548667775,0.468353856,0.495793694,0.450233883,0.435563393,0.403205311,0.427295594,0.359591388,-0.582969908,0.335402222,0.309952966,0.300129552,0.293129426,0.293583798,0.239513564,0.255303141,0.220204387,0.225505408,0.227246057,0.154325388,0.199752451,0.169280343,0.176321289,0.167477769,0.145899825,0.125536605,0.12135882,0.114396878,0.11778856,0.118246595,0.102681041,0.085488536,0.070239495,0.093348714,0.088575438,0.077767637,0.078043344,0.071393113,0.041080937,0.058835859,0.048843231,0.051745413,0.06697217]


# SX-GZO, KRAMER CMAS, 1300C, 10 MINUTES DATA FOR GD
xdata = [12.40596001,17.79057841,23.87262148,30.19653842,36.52045536,42.60249843,47.98711683,54.06925989,60.39307684,66.4752199,71.85973831,78.18365525,84.26569832,90.58961526,96.67175832,102.0559768,108.3802937,114.4626367,120.7859537,126.8682967,132.2527152,138.5770321,144.6583752,150.9826921,156.0811414,162.4054583,168.4868015,174.8111184,181.1344354,186.2328847,192.5572016,198.6395446,204.9628616,211.2871785,216.3856278,222.7089448,228.7912879,235.1156048,240.5000232,246.5823662,252.9056832,258.9880263,265.3123432,270.4107925,276.7341095,283.0584264,289.1397695,295.4640864,300.5625357,306.8868526,313.2101696,319.2925127,325.6168296,330.7152789,337.0385959,343.3629128,349.4442559,354.8296743,360.9110174,367.2353343,373.3176773,379.6409943,385.0254128,391.1077558,397.4320727,403.5134158,409.8377327,415.2221512,421.3044942,427.6278112,433.7101542,439.0945727,445.1769157,451.5012326,457.8245496]
xdata[:] = [x /1000 for x in xdata]
ydata = [9.89231,9.41396,9.2277,8.43706,7.90964,7.63876,7.16472,6.72567,6.23251,5.73609,5.42267,5.03725,4.7584,4.31077,4.14803,3.79136,3.60172,3.26408,2.97397,2.83087,2.56091,2.46369,2.2772,2.09434,1.94178,1.78044,1.57988,1.50183,1.44691,1.33493,1.20359,1.15463,1.0815,0.986342,0.897575,0.736568,0.765653,0.684753,0.663161,0.572472,0.606262,0.479886,0.494207,0.465259,0.438759,0.345015,0.426095,0.400512,0.301205,0.323849,0.248841,0.232741,0.289967,0.190217,0.215721,0.17031,0.244019,0.16404,0.15821,0.209806,0.199097,0.187551,0.137852,0.139485,0.094311,0.104611,0.08307,0.124317,0.194277,0.072871,0.084027,0.093405,0.154589,0.100182,0.143054]

best_fit_data,fit_errors = curve_fit(func, xdata, ydata,[-0.19,10.0,13E-6])

print(best_fit_data)
print(fit_errors)
print()
perr = numpy.sqrt(numpy.diag(fit_errors))
print(perr)

a1,a2,a3 = best_fit_data
print("The best fit parameter values are: ")
print("Co = {}".format(a2))
print("a = {}".format(a1))
print("Do = {}".format(a3))

# Co = 2.86257134475
#a = -0.199139903568
#Do = 6.00843889285e-06

