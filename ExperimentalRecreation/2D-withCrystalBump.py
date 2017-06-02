"""The goal of this code is to do a 2D experimental recreation, similar to GdDiffusionInKramer1300.py
but now with a realistically sized bump in the melt region that represents the physical space taken up
by the crystal. This will include the melt overhangs on both sides showing the incomplete crystal 
coverage of the bottom of the crucible"""

from __future__ import print_function
from fenics import *
from mshr import *
import numpy as np
from subprocess import call

T = 3600			# Final Time
num_steps = 100		# number of time steps
dt = T / num_steps	# time step size
diff = 18E-6		# Diffusivity (mm^2/s)

# Create mesh and define function space
nx = 192
ny = 192
#mesh = RectangleMesh(Point(-2,0), Point(2,4), nx, ny)
#V = FunctionSpace(mesh, 'P', 1)
meltpool = Rectangle(Point(-2,0), Point(2,4))
crystal = Rectangle(Point(-1.4,0), Point(1.4,1))
domain = meltpool - crystal
mesh = generate_mesh(domain, 128)
V = FunctionSpace(mesh, 'P', 1)

# Define the tolerance for the boundary condition
tol = 1.0E-14

# Define the boundary condition
# Including the mixed Neumann and Dirichlet conditions
u_D = Expression('1.0 * (1 - exp(-0.4 * (t/60)))', degree=2, t=0) # POSSIBLE ISSUE IN DEGREE

def boundary_D(x, on_boundary):
	if on_boundary:
		if near(x[0], 0, 1.40) and near(x[1], 1, tol):
			return True
		elif near(x[0], -1.4, tol) and near(x[1], 0.5, 0.5):
			return True
		elif near(x[0], 1.4, tol) and near(x[1], 0.5, 0.5):
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
f = Constant(0)

#F = u*v*dx + dt*diff*dot(grad(u), grad(v))*dx - (u_n + dt*f)*v*dx
a = (u*v + dt*diff*dot(grad(u), grad(v)))*dx
L = (u_n + dt*f)*v*dx

# Create VTK file for saving solution
#call(["ls"])
vtkfile = File('data2D-withCrystalBump/solution.pvd')


# Time-stepping
u = Function(V)
t = 0.0


for n in range(num_steps):

	# Update current time
	t += dt
	u_D.t = t

	# Compute Solution
	solve(a == L, u, bc)

	# Save to file and plot solution
	vtkfile << (u,t)

	# Compute Error at vertices
	# Not done here

	# Print time
	print('t = %.2f' % t)

	# Update previous solution
	u_n.assign(u)

np.savetxt('nodes.csv',mesh.coordinates())