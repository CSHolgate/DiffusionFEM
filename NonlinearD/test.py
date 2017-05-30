"""
FEniCS tutorial demo program: Nonlinear Poisson equation.

  -div(q(u)*grad(u)) = f   in the unit square.
                   u = u_D on the boundary.
"""

from __future__ import print_function

# Warning: from fenics import * will import both `sym` and
# `q` from FEniCS. We therefore import FEniCS first and then
# overwrite these objects.
from fenics import *

def q(u):
	"Return nonlinear coefficient"
	#return 5.0E-6 * (1.0 - 0.36*u + ((0.36*u)**2.0)/2.0)
	return 5.0E-6# * fenics.exp(-0.36*u)

u_D = Expression('3.0 * (1 - exp(-5000 * t))', degree=2, t=0) # POSSIBLE ISSUE IN DEGREE

tol = 1e-14
def boundary(x, on_boundary):
	if on_boundary:
		if near(x[1], 0, tol) and near(x[0],0,0.8):
			return True
		else:
			return False
	else:
		return False

# Create mesh and define function space
mesh = UnitSquareMesh(8, 8)
V = FunctionSpace(mesh, 'P', 1)

bc = DirichletBC(V, u_D, boundary)

# Define variational problem
u = Function(V)  # Note: not TrialFunction!
v = TestFunction(V)
f = Constant(0.0)
F = q(u)*dot(grad(u), grad(v))*dx - f*v*dx

# Compute solution
solve(F == 0, u, bc)

# Save solution to file in VTK format
vtkfile = File('NLPoisson/solution.pvd')
vtkfile << u

# Plot solution
plot(u)

# Compute maximum error at vertices. This computation illustrates
# an alternative to using compute_vertex_values as in poisson.py.
u_e = interpolate(u_D, V)
import numpy as np
error_max = np.abs(u_e.vector().array() - u.vector().array()).max()
print('error_max = ', error_max)

# Hold plot
interactive()
