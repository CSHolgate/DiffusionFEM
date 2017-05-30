"""
Diffusion problem meant to recreate experiments. Has a combination of Dir and Neu boundary
Conditions, just as actual experiments do.
"""

from __future__ import print_function
from fenics import *
import numpy as np
from subprocess import call

T = 3600.0			# Final Time
num_steps = 400		# number of time steps
dt = T / num_steps	# time step size
diff = 5E-6		# Diffusivity (um^2/s)

# Create mesh and define function space
nx = 96
ny = 96
mesh = RectangleMesh(Point(-1,0), Point(1,2), nx, ny)
V = FunctionSpace(mesh, 'P', 1)

# Define the tolerance for the boundary condition
tol = 1.0E-14

# Define the boundary condition
# Including the mixed Neumann and Dirichlet conditions
u_D = Expression('2.65 * (1 - exp(-0.49 * (t/60)))', degree=2, t=0) # POSSIBLE ISSUE IN DEGREE

def boundary_D(x, on_boundary):
	if on_boundary:
		if near(x[0], 0, 0.70) and near(x[1], 0, tol):
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
vtkfile = File('dataDirandNeu/solution.pvd')


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