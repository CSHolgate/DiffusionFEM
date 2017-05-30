"""The goal of this code is to do the time stepping with a while loop that may enable the 
time stepping to stop once a critical concentration is reached in the body center. I've
done something like this in Mathematica and this is just me trying to get it to work here 
well. The output of this should be the time that it took for the system to reach that 
critical value of concentration. I'll start with non-linear 2D with a simple square"""


from __future__ import print_function
from fenics import *
import numpy as np

# Ask the user for information on the diffusivity, and time step
print('What is the desired:')
dt = float(input('Time step? '))
diffarray = [10]#[0.001,0.00316,0.01,0.0316,0.1, 0.316, 1, 3.16, 10, 31.6, 100]#input('Diffusivity? ENTER AS LIST  ')
tfinal = [] # Initiate a list to store times in
afactor = [0.48]#[0.48, 30000] # What exponential factors to calculate

for i in range(0,len(afactor)):
	tfinal.append("a = {}".format(afactor[i]))

	for n in range(0,len(diffarray)):
		#dt = 0.1 		# Time between each time step
		diff = float(diffarray[n])		# Diffusivity (um^2/s)

		# Create mesh and define function space
		#nx = ny = 12
		#mesh = UnitSquareMesh(nx,ny)
		nx = 48
		ny = 48
		mesh = RectangleMesh(Point(-0.5,-0.5), Point(0.5,0.5), nx, ny)
		V = FunctionSpace(mesh, 'P', 1)

		# Define the boundary condition
		u_D = Expression('2.65 * (1 - exp(-a * (t/60)))', degree=2, t=0, a=afactor[i]) # POSSIBLE ISSUE IN DEGREE
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
		vtkfile = File('dataKillAtConc/a{}_d{}/solution.pvd'.format(afactor[i],diffarray[n]))


		# Time-stepping
		u = Function(V)
		t = 0.0

		while u(0,0) < 2.65*0.95:

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
			print('t = %.4f' % t)
			print(u(0,0))

			# Update previous solution
			u_n.assign(u)
			
		# Record time for calculation
		tfinal.append(t)

print(tfinal)
