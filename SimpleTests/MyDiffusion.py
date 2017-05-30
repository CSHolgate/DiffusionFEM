"""
Fenics trial run using some of my diffusion data. 
A test to see how well I am understanding the problem.
"""

from __future__ import print_function
from fenics import *
import numpy as np

T = 2.0				# Final Time
num_steps = 25		# number of time steps
dt = T / num_steps	# time step size
diff = 5E-6		# Diffusivity (um^2/s)

# Create mesh and define function space
#nx = ny = 12
#mesh = UnitSquareMesh(nx,ny)
nx = 48
ny = 96
mesh = RectangleMesh(Point(-1,0), Point(1,4), nx, ny)
V = FunctionSpace(mesh, 'P', 1)

# Define the boundary condition
u_D = Expression('3.0 * (1 - exp(-0.5 * t))', degree=2, t=0) # POSSIBLE ISSUE IN DEGREE
#u_D = Constant(3.0)

def boundary(x, on_boundary):
	return on_boundary

bc = DirichletBC(V, u_D, boundary)

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
vtkfile = File('dataMyDiffusion/solution.pvd')


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
	#print(u(0,0))

	# Update previous solution
	u_n.assign(u)

print(mesh.coordinates())
print(type(mesh.coordinates()))
np.savetxt('nodes.csv',mesh.coordinates(),delimiter=',')
np.savetxt('data.csv',u.vector().array()[::-1],delimiter='\t')
