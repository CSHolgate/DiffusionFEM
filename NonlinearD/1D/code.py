"""
Fenics program on an attempt to create a way to model concentration dependent
diffusivity. This will be a NON-LINEAR solver
"""

from __future__ import print_function
from fenics import *
import fenics
import numpy
from subprocess import call

T = 3480.0
num_steps = 300
dt = T/num_steps

# Define the non-linear coefficient Diff
def Diff(u):
	"Return nonlinear coefficient"
	#return 5.0E-6 * (1.0 - 0.36*u + ((0.36*u)**2.0)/2.0)
	return 6.215629E-6 * fenics.exp(-0.19914*u)



# Create mesh and define function space
mesh = IntervalMesh(10000, 0, 2)
V = FunctionSpace(mesh, 'P', 1)

# Define the tolerance of the boundary condition
tol = 1.0E-14

# Define the boundary condition
# Including mixed Neumann (0 at boundaries) and Dirichlet conditions
u_D = Expression('2.8625716 * (1 - exp(-30000 * t))', degree=2, t=0) # POSSIBLE ISSUE IN DEGREE

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

#F = (u*v + dt*1.0*dot(grad(u), grad(v)))*dx - (u_n + dt*f)*v*dx
#a = (u*v + dt*diff*dot(grad(u), grad(v)))*dx
#L = (u_n + dt*f)*v*dx

# Create VTK file for saving solution
vtkfile = File('dataDiffFConc/solution.pvd')

# Time stepping
#u = Function(V)
t = 0.0

for n in range(num_steps):

	# Update current time
	t += dt
	u_D.t = t

	# Compute solution
	solve(F == 0, u, bc)

	# Save to file and plot solution
	vtkfile << (u,t)

	# Print time
	print('t= %2.f' % t)

	# Update previous solution
	u_n.assign(u)

#print(u.vector().array())
