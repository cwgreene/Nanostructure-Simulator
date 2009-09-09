"""This demo program solves Poisson's equation

    - div grad u(x, y) = f(x, y)

on the unit square with source f given by

    f(x, y) = 500*exp(-((x - 0.5)^2 + (y - 0.5)^2) / 0.02)

and boundary conditions given by

    u(x, y) = 0 for x = 0 or x = 1
"""

__author__ = "Anders Logg (logg@simula.no)/Chris Greene"
__date__ = "2007-08-16 -- 2008-12-13 / 2009-08-18"
__copyright__ = "Copyright (C) 2007-2008 Anders Logg"
__license__  = "GNU LGPL Version 2.1"

from dolfin import *
import montecarlo_mockup as mc
import numpy as np
import dolfin_util as du
import time
#dolfin_set("linear algebra backend","Epetra"

def custom_func(mesh,V,particles):
	f = Function(V)
	for p in particles:
		du.alter_cellid(mesh,f,p.id,p.charge)
	return f

# Create mesh and define function space
meshSizeX = 50
meshSizeY = 50
mesh = UnitSquare(meshSizeX,meshSizeY)
#plot(mesh)
V = FunctionSpace(mesh, "Lagrange", 2)

# Define Dirichlet boundary (x = 0 or x = 1)
class DirichletBoundary(SubDomain):
    def inside(self, x, on_boundary):
        return (x[0] < DOLFIN_EPS  or (x[0] > 1.0 - DOLFIN_EPS) or 
		x[1] < DOLFIN_EPS or (x[1] > 1.0 -DOLFIN_EPS))

# Define boundary condition
u0 = Constant(mesh, 0.0)
bc = DirichletBC(V, u0, DirichletBoundary())

# Define variational problem
v = TestFunction(V)
u = TrialFunction(V)

#function
particles = []
for x in np.arange(0,1.,1./meshSizeX):
	for y in np.arange(0.,1.,1./meshSizeY):
		#electrons
		particles += mc.init_electrons(1,[[x,y]],charge=-10,mesh=mesh)
		#holes
		particles += mc.init_electrons(1,[[x,y]],charge=10,mesh=mesh)
f = custom_func(mesh,V,particles)
file = File("poisson_attract.pvd")
dfile = File("density_attract.pvd")
for x in range(200):
	u = TrialFunction(V)
	v = TestFunction(V)

	a = dot(grad(v), grad(u))*dx
	L = v*f*dx
	# Compute solution
	problem = VariationalProblem(a, L, bc)
	sol = problem.solve()
	# Plot solution
	file << sol
	dfile << f
	print "Starting Step ",x
	start = time.time()
	mc.MonteCarlo(mesh,sol,f,particles)
	print "Took: ",time.time()-start
#	plot(f)
	du.delete(problem)
file << sol
dfile << f
# Hold plot
interactive()
